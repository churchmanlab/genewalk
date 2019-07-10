import re
import os
import pickle
import logging
import itertools
import pandas as pd
import networkx as nx
from indra.databases import go_client
from goatools.obo_parser import GODag
from genewalk.resources import get_go_obo, get_goa_gaf, get_pc
from genewalk.get_indra_stmts import get_famplex_links_from_stmts

logger = logging.getLogger('genewalk.nx_mg_assembler')


# TODO: these assemblers have a lot of code duplications, they could be
# refactored to derive from a single assembler class

def load_network(network_type, network_file, genes):
    if network_type == 'pc':
        MG = PcNxMgAssembler(genes)
    elif network_type == 'indra':
        logger.info('Loading %s' % network_file)
        with open(network_file, 'rb') as fh:
            stmts = pickle.load(fh)
        MG = IndraNxMgAssembler(stmts)
    elif network_type == 'edge_list':
        logger.info('Loading user-provided GeneWalk Network from %s.' %
                    network_file)
        MG = UserNxMgAssembler(network_file, gwn_format='el')
    elif network_type == 'sif':
        logger.info('Loading user-provided GeneWalk Network from %s.' %
                    network_file)
        MG = UserNxMgAssembler(network_file, gwn_format='sif')
    else:
        raise ValueError('Unknown network_type: %s' % network_type)
    return MG


def _load_goa_gaf():
    goa_ec = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP',
              'HGI', 'HEP', 'IBA', 'IBD'}
    goa = pd.read_csv(get_goa_gaf(), sep='\t', skiprows=23, dtype=str,
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


def _build_go_ontology():
    go_ontology = {}
    for go_term in go_dag.values():
        if go_term.is_obsolete:
            continue
        for parent_term in go_term.parents:
            if parent_term.is_obsolete:
                continue
            if go_term.id in go_ontology:
                go_ontology[go_term.id].append(parent_term.id)
            else:
                go_ontology[go_term.id] = [parent_term.id]
    return go_ontology


go_dag = GODag(get_go_obo())
go_ontology = _build_go_ontology()
goa = _load_goa_gaf()


class NxMgAssembler(object):
    def __init__(self, genes):
        self.genes = genes
        self.graph = nx.MultiGraph()

    @staticmethod
    def _get_go_terms_for_gene(gene):
        # Filter to rows with the given gene's UniProt ID
        df = goa[goa['DB_ID'] == gene['UP']]
        go_ids = sorted(list(set(df['GO_ID'])))
        return go_ids

    def add_go_annotations(self):
        logger.info('Adding GO annotations for genes to graph.')
        for gene in self.genes:
            go_ids = self._get_go_terms_for_gene(gene)
            for go_id in go_ids:
                go_term = go_dag[go_id]
                if go_term.is_obsolete:
                    continue
                self.graph.add_node(go_term.id,
                                    name=go_term.name.replace(' ', '_'),
                                    GO=go_term.id, source='go')
                # TODO: do we need qualifiers here as labels?
                self.graph.add_edge(gene['HGNC_SYMBOL'], go_term.id,
                                    label='assoc_with')

    def add_go_ontology(self):
        """Add to self.graph the GO ontology (GO:IDs and their relations) in
        the form of labeled edge (relation type, eg is_a) and new nodes
        (GO:IDs).
        """
        logger.info('Adding GO ontology edges to graph.')
        for go_term in go_dag.values():
            if go_term.is_obsolete:
                continue
            self.graph.add_node(go_term.id,
                                name=go_term.name.replace(' ', '_'),
                                GO=go_term.id, source='go')
            for parent_term in go_term.parents:
                if parent_term.is_obsolete:
                    continue
                self.graph.add_node(go_term.id,
                                    name=go_term.name.replace(' ', '_'),
                                    GO=go_term.id, source='go')
                self.graph.add_edge(go_term.id, parent_term.id,
                                    label='GO:is_a')

    def node2edges(self, node_key):
        return self.graph.edges(node_key, keys=True)

    def save_graph(self, fname):
        nx.write_graphml(self.graph, fname)


class PcNxMgAssembler(NxMgAssembler):
    """The PcNxMgAssembler assembles a GeneWalk Network with gene reactions
    from Pathway Commons and GO ontology and annotations into a networkx
    (undirected)  MultiGraph including edge attributes.

    Parameters
    ----------
    genes : list

    Attributes
    ----------
    graph : networkx.MultiGraph
        A GeneWalk Network that is assembled by this assembler.
    GOA : pandas.DataFrame
        GO annotation in pd.dataframe format
    OGO : goatools.GODag
        GO ontology, GODag object (see goatools) 
    """
    def __init__(self, genes):
        super().__init__(genes)
        self.add_go_ontology()
        self.add_go_annotations()
        self.add_pc_edges()

    def add_pc_edges(self):
        """Assemble a nx.MultiGraph from the Pathway Commons sif file
        (nodeA <relationship type> nodeB).
        """
        logger.info('Adding gene edges from Pathway Commons to graph.')
        gwn_df = pd.read_csv(get_pc(), sep='\t', dtype=str, header=None)
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
        pc_sub = pc.subgraph(hgnc_symbols)
        gene2hgnc_dict = dict(zip(hgnc_symbols, hgnc_ids))
        nx.set_node_attributes(pc_sub, gene2hgnc_dict, 'HGNC')
        gene2up_dict = dict(zip(hgnc_symbols, up_ids))
        nx.set_node_attributes(pc_sub, gene2up_dict, 'UP')
        # make a copy to unfreeze graph
        pc_graph = nx.MultiGraph(pc_sub)
        logger.info('Number of PC originating nodes %d' %
                    nx.number_of_nodes(pc_graph))
        self.graph = nx.compose(self.graph, pc_graph)


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
    def __init__(self, genes, stmts):
        self.stmts = stmts
        super().__init__(genes)
        self.add_indra_edges()
        self.add_fplx_edges()
        self.indra_nodes = {}

    def add_indra_edges(self):
        """Assemble the graph from the assembler's list of INDRA Statements. 
        Edge attribute are given by statement type and index in list of stmts
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
            edge_key = '%d_%s' % (i, type(st).__name__)
            # Iterate over all the agent combinations and add edge
            for a, b in itertools.combinations(agents, 2):
                a_node = self.add_agent_node(a)
                b_node = self.add_agent_node(b)
                self.graph.add_edge(a_node, b_node, label=edge_key)
        logger.info('Number of INDRA originating nodes %d.' %
                    len(self.indra_nodes))


    def add_fplx_edges(self):
        """Add to self.graph an edge (label: 'FPLX:is_a') between the gene
        family to member annotation edges.
        """
        links = get_famplex_links_from_stmts(self.stmts)
        for s, t in links:
            self.graph.add_edge(s, t, label='FPLX:is_a')

    def add_agent_node(self, agent):
        go_id = agent.db_refs.get('GO')
        if go_id:
            go_id = go_id if go_id.startswith('GO:') else 'GO:%s' % go_id
            node_key = go_id
            name = go_client.get_go_label(go_id)
            self.graph.add_node(node_key, name=name.replace(' ', '_'),
                                source='indra', **agent.db_refs)
        else:
            node_key = agent.name
            self.graph.add_node(node_key, name=agent.name, **agent.db_refs,
                                source='indra')
        self.indra_nodes.add(node_key)
        return node_key

    def node2stmts(self, node_key):
        matching_stmts = []
        node_name = self.graph.node[node_key]['name']
        for stmt in self.stmts:
            for agent in stmt.agent_list():
                if agent is not None:
                    agent_name = agent.name
                    if agent_name == node_name:
                        matching_stmts.append(stmt)
                        break
        return matching_stmts


class UserNxMgAssembler(object):
    """The UserNxMgAssembler loads a user-provided GeneWalk Network from
    file.

    Parameters
    ----------
    filepath : Optional[str]
        Path to the user-provided genewalk network file, assumed to contain
        gene symbols and GO:IDs. See gwn_format for supported format details.
    gwn_format : Optional[str]
        'el' (default, edge list: nodeA nodeB (if more columns
        present: interpreted as edge attributes) \
        or 'sif' (simple interaction format: nodeA <relationship type> nodeB).
        Do not include column headers.

    Attributes
    ----------
    graph : networkx.MultiGraph
        A GeneWalk Network that is loaded by this assembler.
    """
    # TODO: reimplement this as a method in main assembler class?
    # TODO: gwn_format is unused here
    def __init__(self, filepath='~/genewalk/gwn.txt', gwn_format='el'):
        self.graph = nx.MultiGraph()
        self.filepath = filepath
        
    def MG_from_file(self):
        """Assemble the GeneWalk Network from the user-provided file path."""
        gwn_df = pd.read_csv(self.filepath, dtype=str, header=None)
        col_mapper = {}
        if self.gwn_format == 'el':
            col_mapper[0] = 'source'
            col_mapper[1] = 'target'
            if len(gwn_df.columns) > 2:
                col_mapper[2] = 'rel_type'
                if len(gwn_df.columns) > 3:
                    for c in gwn_df.columns[3:]:
                        col_mapper[c] = 'edge_attr'+str(c-1)
                edge_attributes = True
            else:
                edge_attributes = False
        elif self.gwn_format == 'sif':
            col_mapper[0] = 'source'
            col_mapper[1] = 'rel_type'
            col_mapper[2] = 'target'
            edge_attributes = True
            
        gwn_df.rename(mapper=col_mapper,axis='columns')
        self.graph = nx.from_pandas_edgelist(gwn_df, 'source', 'target',
                                             edge_attr=edge_attributes,
                                             create_using=nx.MultiGraph)

    def node2edges(self, node_key):
        return self.graph.edges(node_key,keys=True)        
