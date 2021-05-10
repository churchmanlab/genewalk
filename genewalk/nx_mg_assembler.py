import pickle
import logging
import itertools
import pandas as pd
import networkx as nx
from goatools.obo_parser import GODag
from genewalk.resources import ResourceManager

logger = logging.getLogger('genewalk.nx_mg_assembler')


def load_network(network_type, network_file, genes, resource_manager=None):
    """Return a network assembler of the given type based on a set of genes.

    Parameters
    ----------
    network_type : str
        The type of the network to be constructed.
    network_file : str
        The path to a file containing information to construct the network.
    genes : list
        A list of gene references.
    resource_manager : Optional[:py:class:`genewalk.resources.ResourceManager`]
        A resource manager object which, if specified, is used to get the
        resource files. Otherwise, the default resource manager is used.

    Returns
    -------
    :py:class:`genewalk.nx_mg_assembler.NxMgAssembler`
        An instance of an NxMgAssembler containing the assembled networkx
        MultiGraph as its graph attribute.
    """
    if not resource_manager:
        resource_manager = None
    if network_type == 'pc':
        mg = PcNxMgAssembler(genes, resource_manager=resource_manager)
    elif network_type == 'indra':
        logger.info('Loading %s' % network_file)
        with open(network_file, 'rb') as fh:
            stmts = pickle.load(fh)
        mg = IndraNxMgAssembler(genes, stmts,
                                resource_manager=resource_manager)
    elif network_type in {'edge_list', 'sif', 'sif_annot', 'sif_full'}:
        logger.info('Loading user-provided GeneWalk Network from %s.' %
                    network_file)
        mg = UserNxMgAssembler(genes, resource_manager,
                               network_file, gwn_format=network_type)
    else:
        raise ValueError('Unknown network_type: %s' % network_type)
    return mg


class NxMgAssembler(object):
    """Class which assembles a networkx MultiGraph based on a list of genes.

    Parameters
    ----------
    genes : list of dict
        A list of gene references based on which the graph is assembled.

    Attributes
    ----------
    graph : networkx.MultiGraph
        The assembled graph containing links for interactions between genes,
        GO annotations for genes, and the GO ontology.
    """

    def __init__(self, genes, resource_manager=None):
        self.genes = genes
        self.graph = nx.MultiGraph()
        if not resource_manager:
            self.resource_manager = ResourceManager()
        else:
            self.resource_manager = resource_manager
        self.go_dag = None
        self.goa = None

    def get_go_annots(self):
        if self.goa is None:
            self.goa = self._load_goa_gaf()
        return self.goa

    def get_go_dag(self):
        if self.go_dag is None:
            self.go_dag = GODag(self.resource_manager.get_go_obo())
        return self.go_dag

    def _get_go_terms_for_gene(self, gene):
        # Filter to rows with the given gene's UniProt ID
        if ('UP' not in gene) or ('HGNC_SYMBOL' not in gene):
            return []
        elif gene['HGNC_SYMBOL'] not in self.graph:
            return []
        goa = self.get_go_annots()
        df = goa[goa['DB_ID'] == gene['UP']]
        go_ids = sorted(list(set(df['GO_ID'])))
        return go_ids

    def add_go_annotations(self):
        """Add edges between gene nodes and GO nodes based on GO
        annotations."""
        logger.info('Adding GO annotations for genes in graph.')
        go_dag = self.get_go_dag()
        for gene in self.genes:
            go_ids = self._get_go_terms_for_gene(gene)
            for go_id in go_ids:
                if go_id in go_dag:
                    go_term = go_dag[go_id]
                    if go_term.is_obsolete:
                        continue
                    self.graph.add_node(go_term.id,
                                        name=go_term.name,
                                        GO=go_term.id,
                                        domain=go_term.namespace)
                    self.graph.add_edge(gene['HGNC_SYMBOL'], go_term.id,
                                        label='GO:annotation')

    def add_go_ontology(self):
        """Add edges between GO nodes based on the GO ontology."""
        logger.info('Adding GO ontology edges to graph.')
        go_dag = self.get_go_dag()
        for go_term in list(go_dag.values()):
            if go_term.is_obsolete:
                continue
            self.graph.add_node(go_term.id,
                                name=go_term.name,
                                GO=go_term.id,
                                domain=go_term.namespace)
            for parent_term in go_term.parents:
                if parent_term.is_obsolete:
                    continue
                self.graph.add_node(go_term.id,
                                    name=go_term.name,
                                    GO=go_term.id,
                                    domain=go_term.namespace)
                self.graph.add_edge(go_term.id, parent_term.id,
                                    label='GO:is_a')

    def node2edges(self, node_key):
        """Return the edges corresponding to a node."""
        return self.graph.edges(node_key, keys=True)

    def save_graph(self, fname):
        """Save the file into a GraphML file.

        Parameters
        ----------
        fname : str
            The name of the file to save the graph into.
        """
        nx.write_graphml(self.graph, fname)

    def _load_goa_gaf(self):
        """Load the gene/GO annotations as a pandas data frame."""
        goa_ec = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA',
                  'HMP', 'HGI', 'HEP', 'IBA', 'IBD'}
        goa = pd.read_csv(self.resource_manager.get_goa_gaf(), sep='\t',
                          comment='!', dtype=str,
                          header=None,
                          names=['DB',
                                 'DB_ID',
                                 'DB_Symbol',
                                 'Qualifier',
                                 'GO_ID',
                                 'DB_Reference',
                                 'Evidence_Code',
                                 'With_From',
                                 'Aspect',
                                 'DB_Object_Name',
                                 'DB_Object_Synonym',
                                 'DB_Object_Type',
                                 'Taxon',
                                 'Date',
                                 'Assigned',
                                 'Annotation_Extension',
                                 'Gene_Product_Form_ID'])
        goa = goa.sort_values(by=['DB_ID', 'GO_ID'])
        # Filter out all "NOT" negative evidences
        goa['Qualifier'].fillna('', inplace=True)
        goa = goa[~goa['Qualifier'].str.startswith('NOT')]
        # Filter to rows with evidence code corresponding to experimental
        # evidence
        goa = goa[goa['Evidence_Code'].isin(goa_ec)]
        return goa


class PcNxMgAssembler(NxMgAssembler):
    """The PcNxMgAssembler assembles a GeneWalk Network with gene reactions
    from Pathway Commons and GO ontology and annotations into a networkx
    (undirected)  MultiGraph including edge attributes.

    Parameters
    ----------
    genes : list fo dict
        A list of gene references based on which the network is assembled.

    Attributes
    ----------
    graph : networkx.MultiGraph
        A GeneWalk Network that is assembled by this assembler.
    """
    def __init__(self, genes, resource_manager=None):
        super().__init__(genes, resource_manager)
        self.add_pc_edges()
        self.add_go_annotations()
        self.add_go_ontology()

    def add_pc_edges(self):
        """Add edges between gene nodes based on PathwayCommons
        interactions."""
        logger.info('Adding gene edges from Pathway Commons to graph.')
        gwn_df = pd.read_csv(self.resource_manager.get_pc(), sep='\t',
                             dtype=str, header=None)
        col_mapper = {}
        col_mapper[0] = 'source'
        col_mapper[1] = 'rel_type'
        col_mapper[2] = 'target'
        edge_attributes = True
        gwn_df = gwn_df.rename(mapper=col_mapper, axis='columns')
        pc = nx.from_pandas_edgelist(gwn_df, source='source', target='target',
                                     edge_attr=edge_attributes,
                                     create_using=nx.MultiGraph)
        # subset over genes in the input gene list
        hgnc_symbols = [g['HGNC_SYMBOL'] for g in self.genes]
        hgnc_ids = [g['HGNC'] for g in self.genes]
        up_ids = [g['UP'] for g in self.genes]
        pc.add_nodes_from(hgnc_symbols)
        pc_sub = pc.subgraph(hgnc_symbols)
        gene2hgnc_dict = dict(zip(hgnc_symbols, hgnc_ids))
        nx.set_node_attributes(pc_sub, gene2hgnc_dict, 'HGNC')
        gene2up_dict = dict(zip(hgnc_symbols, up_ids))
        nx.set_node_attributes(pc_sub, gene2up_dict, 'UP')
        # Make a copy to unfreeze graph
        self.graph = nx.MultiGraph(pc_sub)
        logger.info('Number of PC originating nodes %d' %
                    nx.number_of_nodes(self.graph))


class IndraNxMgAssembler(NxMgAssembler):
    """The IndraNxMgAssembler assembles INDRA Statements and GO ontology /
    annotations into a networkx (undirected) MultiGraph including edge
    attributes. This code is based on INDRA's SifAssembler
    http://indra.readthedocs.io/en/latest/_modules/indra/assemblers/sif_assembler.html

    Parameters
    ----------
    stmts : list[indra.statements.Statement]
        A list of INDRA Statements to be added to the assembler's list
        of Statements.

    Attributes
    ----------
    graph : networkx.MultiGraph
        A GeneWalk Network that is assembled by this assembler.
    """
    def __init__(self, genes, stmts, resource_manager=None):
        self.indra_nodes = set()
        self.stmts = stmts
        super().__init__(genes, resource_manager)
        self.add_indra_edges()
        self.add_fplx_edges()
        self.add_go_annotations()
        self.add_go_ontology()

    def add_indra_edges(self):
        """Add edges between gene nodes and GO nodes based on INDRA Statements.
        """
        logger.info('Adding nodes from INDRA statements.')
        for i, st in enumerate(self.stmts):
            # Get all agents in the statement
            agents = [a for a in st.agent_list() if a is not None]
            # Only include edges for statements with at least 2 Agents
            # excludes (irrelevant) stmt types: Translocation, ActiveForm,
            # SelfModification
            if len(agents) < 2:
                continue
            # Create a label that is unique to the statement and its type
            edge_type = type(st).__name__
            edge_key = '%d_%s' % (i, edge_type)
            # Iterate over all the agent combinations and add edge
            for a, b in itertools.combinations(agents, 2):
                a_node = self.add_agent_node(a)
                b_node = self.add_agent_node(b)
                self.graph.add_edge(a_node, b_node, key=edge_key,
                                    label=edge_type)
        hgnc_symbols = [g['HGNC_SYMBOL'] for g in self.genes]
        hgnc_ids = [g['HGNC'] for g in self.genes]
        up_ids = [g['UP'] for g in self.genes]
        self.graph.add_nodes_from(hgnc_symbols)
        gene2hgnc_dict = dict(zip(hgnc_symbols, hgnc_ids))
        nx.set_node_attributes(self.graph, gene2hgnc_dict, 'HGNC')
        gene2up_dict = dict(zip(hgnc_symbols, up_ids))
        nx.set_node_attributes(self.graph, gene2up_dict, 'UP')
        logger.info('Number of INDRA originating nodes %d.' %
                    len(self.indra_nodes))

    def add_fplx_edges(self):
        """Add edges between gene nodes and families/complexes they are part
        of."""
        from genewalk.get_indra_stmts import get_famplex_links_from_stmts
        links = get_famplex_links_from_stmts(self.stmts)
        for s, t in links:
            self.graph.add_edge(s, t, label='FPLX:is_a')

    def add_agent_node(self, agent):
        """Add a node corresponding to an INDRA Agent."""
        go_id = agent.db_refs.get('GO')
        go_dag = self.get_go_dag()
        if go_id:
            go_id = go_id if go_id.startswith('GO:') else 'GO:%s' % go_id
            node_key = go_id
            # INDRA standardizes GO names so this is generally not
            # necessary
            try:
                name = go_dag[go_id].name
            except KeyError:
                name = agent.name
            self.graph.add_node(node_key, name=name,
                                source='indra', **agent.db_refs)
        else:
            node_key = agent.name
            self.graph.add_node(node_key, name=agent.name, **agent.db_refs,
                                source='indra')
        self.indra_nodes.add(node_key)
        return node_key

    def node2stmts(self, node_key):
        """Return the INDRA Statements given the key of a graph node."""
        matching_stmts = []
        node_name = self.graph.nodes[node_key]['name']
        for stmt in self.stmts:
            for agent in stmt.agent_list():
                if agent is not None:
                    agent_name = agent.name
                    if agent_name == node_name:
                        matching_stmts.append(stmt)
                        break
        return matching_stmts


class UserNxMgAssembler(NxMgAssembler):
    """Loads a user-provided GeneWalk Network from a given file.

    Parameters
    ----------
    filepath : str
        Path to the user-provided genewalk network file, assumed to contain
        gene symbols and GO IDs. See gwn_format for supported format details.
    gwn_format : Optional[str]
        'el' (default, edge list: nodeA nodeB (if more columns
        present: interpreted as edge attributes) \
        or 'sif' (simple interaction format: nodeA,<relationship type>,nodeB).
        Do not include column headers.

    Attributes
    ----------
    graph : networkx.MultiGraph
        A GeneWalk Network that is loaded by this assembler.
    """
    def __init__(self, genes, resource_manager, filepath,
                 gwn_format='el'):
        super().__init__(genes, resource_manager=resource_manager)
        self.graph = nx.MultiGraph()
        self.filepath = filepath
        self.gwn_format = gwn_format
        self.add_network_edges()

    def add_network_edges(self):
        """Assemble the GeneWalk Network from the user-provided file path."""
        gwn_df = pd.read_csv(self.filepath, dtype=str, header=None)
        col_mapper = {}
        if self.gwn_format == 'edge_list':
            col_mapper[0] = 'source'
            col_mapper[1] = 'target'
            edge_attributes = None
        elif self.gwn_format in {'sif', 'sif_annot', 'sif_full'}:
            col_mapper[0] = 'source'
            col_mapper[1] = 'rel_type'
            col_mapper[2] = 'target'
            edge_attributes = True
        else:
            raise ValueError('%s is not a valid GeneWalk network format'
                             % self.gwn_format)

        gwn_df.rename(mapper=col_mapper, axis='columns',
                      inplace=True)
        self.graph = nx.from_pandas_edgelist(gwn_df, 'source', 'target',
                                             edge_attr=edge_attributes,
                                             create_using=nx.MultiGraph)
        logger.info('The graph loaded from %s contains %d nodes'
                    ' including %d GO terms' %
                    (self.filepath, len(self.graph),
                     len([n for n in self.graph if n.startswith('GO:')])))
        gene_list_genes = \
            {(g['ID'] if 'ID' in g else g['HGNC_SYMBOL']) for g in self.genes}
        non_gene_list_non_go_nodes = \
            [n for n in self.graph if not n.startswith('GO:')
             and n not in gene_list_genes]
        if non_gene_list_non_go_nodes:
            logger.info('Removing %d gene nodes from input network '
                        'since they are not in the input gene list.' %
                        len(non_gene_list_non_go_nodes))
        self.graph.remove_nodes_from(non_gene_list_non_go_nodes)
        # If the GO annotations are not provided as part of the SIF
        # then we add those
        if self.gwn_format in {'sif', 'edge_list'}:
            self.add_go_annotations()
        # If the GO DAG is not provided as part of the SIF then we add
        # it
        if self.gwn_format in {'sif', 'sif_annot', 'edge_list'}:
            self.add_go_ontology()
        # If the SIF contains everything then we still have to add
        # some basic node meta-data to the GO nodes for later steps
        if self.gwn_format == 'sif_full':
            go_dag = self.get_go_dag()
            for node in self.graph.nodes:
                if node.startswith('GO:'):
                    go_term = go_dag.get(node)
                    if go_term:
                        self.graph.nodes[node]['GO'] = go_term.id
                        self.graph.nodes[node]['name'] = go_term.name
                        self.graph.nodes[node]['domain'] = go_term.namespace
