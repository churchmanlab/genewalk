import json
from builtins import dict, str
import numpy
import itertools
import networkx as nx
import re
import pandas as pd
from indra.statements import *
from goatools.obo_parser import GODag

class Nx_MG_Assembler(object):
    """The Nx_MG_Assembler assembles INDRA Statements and GO ontology /
    annotations into a networkx (undirected) MultiGraph including edge
    attributes. This code is based on INDRA's SifAssembler
    http://indra.readthedocs.io/en/latest/_modules/indra/assemblers/sif_assembler.html

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be added to the assembler's list
        of Statements.
    GOpath : str
        Path to the goa_human.gaf file

    Attributes
    ----------
    graph : networkx.MultiGraph
        A networkx graph that is assembled by this assembler.
    GOA : pandas.DataFrame
        GO annotation in pd.dataframe format
    OGO : goatools.GODag
        GO ontology, GODag object (see goatools)
    """
    def __init__(self,stmts=None,GOpath='~/'):
        self.stmts = [] if stmts is None else stmts
        self.graph = nx.MultiGraph()
        self.GOpath = GOpath
        self.GOA = pd.read_csv(self.GOpath+'goa_human.gaf', sep='\t',skiprows=23,dtype=str,header=None, 
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
        self.GOA = self.GOA.sort_values(by=['DB_ID','GO_ID'])
        self.OGO = GODag(GOpath+'go.obo')#dict
        self.EC_GOA=['EXP','IDA','IPI','IMP','IGI','IEP','HTP','HDA','HMP','HGI','HEP','IBA','IBD']
    
    def MG_from_INDRA(self):
        """Assemble the graph from the assembler's list of INDRA Statements. 
        Edge attribute are given by statement type and index in list of stmts
        """
        N=len(self.stmts)
        for i in range(N):
            if i%1000 == 0:
                print(i,"/",N)
            st = self.stmts[i]
            # Get all agents in the statement
            agents = st.agent_list()
            # Filter out None Agent
            agents = [a for a in agents if a is not None]
            # Only include edges for statements with at least 2 Agents
            if len(agents) < 2:#excludes (irrelevant) stmt types: Translocation, ActiveForm, SelfModification
                continue
            edge_attr = str(i)+'_'+type(st).__name__
            #print(edge_attr)
            #Iterate over all the agent combinations and add edge
            for a, b in itertools.combinations(agents, 2):
                self._add_INnode_edge(a, b, edge_attr)
    
    def add_FPLXannotations(self,filename):
        """Add to self.graph an edge (label: 'FPLX:is_a') between the gene family to member annotation edges.
   
        Parameters
        ----------
        filename : str specifying the .csv file with list of tuples with the first element of the tuple
        a child gene name or FamPlex entry, and the second element
        a parent FamPlex entry, e.g. ('KRAS', 'RAS')
        """
        FPLX=pd.read_csv(filename,sep=',',dtype=str,header=None)
        # Add protein family/complex links
        for i in FPLX.index:
            s=FPLX[0][i]
            t=FPLX[1][i]
            edge_attr='FPLX:is_a'
            self._add_edge(s, t, edge_attr)

    def add_GOannotations(self):
        """Add to self.graph the GO annotations (GO:IDs) of proteins (ie, the
        subset of self.graph nodes that contain UniprotKB:ID) in the form of
        labeled edges (see _GOA_from_UP for details) and new nodes (GO:IDs).
        """
        IN_nodes=list(nx.nodes(self.graph))
        N=len(IN_nodes)
        j=0#counter for duration
        for n in IN_nodes: 
            if j%100 == 0:
                print(j,"/",N)
            if 'UP' in self.graph.node[n].keys():#node is INDRA gene/protein
                UP=self.graph.node[n]['UP']
                GOan=self._GOA_from_UP(UP)
                for i in GOan.index:
                    GOID=GOan['GO_ID'][i]
                    eattr=GOan['Qualifier'][i]
                    if self.OGO[GOID].is_obsolete == False:
                        self._add_GOnode(GOID,'0')
                        self._add_edge(n,GOID,eattr)
            j=j+1
                     
    def add_GOontology(self):
        """Add to self.graph the GO ontology (GO:IDs and their relations) in
        the form of labeled edge (relation type, eg is_a) and new nodes
        (GO:IDs).
        """
        for goid in self.OGO.keys():
            GOT=self.OGO[goid]
            if GOT.is_obsolete == False:
                self._add_GOnode(GOT.id,'0')
                for pa in GOT.parents:
                    if pa.is_obsolete == False:
                        self._add_GOnode(pa.id,'0')
                        self._add_edge(GOT.id,pa.id,'GO:is_a')
    
    def _GOA_from_UP(self,UP):
        SEL=self.GOA[self.GOA['DB_ID']==UP][['GO_ID','Evidence_Code']].drop_duplicates()#UP matching GOIDs and Qualif
        for i in SEL.index:
            if (SEL['Evidence_Code'][i] in self.EC_GOA) & (SEL['GO_ID'][i] in self.OGO):
                pass
            else:
                SEL=SEL.drop(i)#Insufficient evidence for annotation or not present in OGO: obsolete GO:ID, so drop.
        SEL.insert(loc=1,column='Qualifier', value=pd.Series('GOan', index=SEL.index))#add new column
        return SEL.drop(columns=['Evidence_Code'])
 
    def _add_INnode_edge(self, s, t, attributes):
        if s is not None:
            s = self._add_INnode(s)
            t = self._add_INnode(t)
            self._add_edge(s, t, attributes)

    def _add_INnode(self, ag):
        if 'GO' in ag.db_refs:
            node_key=ag.db_refs['GO']
            if re.search(r'GO:', node_key) is None:#double check if GO: is present in GO:ID
                node_key='GO:'+node_key
            self._add_GOnode(node_key,'1')
            for attr in ag.db_refs.keys():#copy over any other identifiers (UP,HGNC,GO,TXT,ChEBI etc) as node attr.
                if attr != 'GO':
                    self.graph.node[node_key][attr]=ag.db_refs[attr]
        else:
            node_key = ag.name
            self.graph.add_node(node_key,name=node_key,INDRA='1')
            for attr in ag.db_refs.keys():#copy over the identifiers (UP,HGNC,TXT,ChEBI etc) as node attribute
                self.graph.node[node_key][attr]=ag.db_refs[attr]
        return node_key

    def _add_GOnode(self,GOID,indra):
        GOT=self.OGO[GOID]
        nameGO=GOT.name 
        nameGO=nameGO.replace(" ","_") 
        self.graph.add_node(GOID,name=nameGO,GO=GOID)#nx ensures no duplicate nodes with same key will be created
        if 'INDRA' not in self.graph.node[GOID].keys():#not yet present, so assign origin: INDRA or GOA/OGO
            self.graph.node[GOID]['INDRA']=indra
        
    def _add_edge(self, s, t, edge_attributes=None):
        if edge_attributes is None:
            self.graph.add_edge(s, t, label='NA')
        else:
            self.graph.add_edge(s, t, label=edge_attributes)
        
    def node2stmts(self, node_key):
        matching_stmts = []
        node_name=self.graph.node[node_key]['name']
        for stmt in self.stmts:
            for agent in stmt.agent_list():
                if agent is not None:
                    agent_name = agent.name
                    if agent_name == node_name:
                        matching_stmts.append(stmt)
                        break
        return matching_stmts
    
    def node2edges(self, node_key):
        return self.graph.edges(node_key,keys=True)
    
    def save_graph(self,folder='~/',filename='test'):
        nx.write_graphml(self.graph,folder+filename+'.xml')
