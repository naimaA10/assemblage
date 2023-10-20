#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "AMMICHE Naïma"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["AMMICHE Naïma"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "AMMICHE Naïma"
__email__ = "naima.ammiche@gmail.com"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    file = open(fastq_file, "r")
    for line in file :
        yield next(file).strip()
        next(file)
        next(file)


def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    for i in range(len(read)-(kmer_size-1)) : 
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    frequence_kmer = {}
    for sequence in read_fastq(fastq_file): 
        for kmer in cut_kmer(sequence, kmer_size) :
            if kmer not in frequence_kmer :
                frequence_kmer[kmer]=1
            else :
                frequence_kmer[kmer]+=1
    return frequence_kmer

def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    graphe = nx.DiGraph()
    for kmer in kmer_dict.keys():
        prefixe = kmer[0:-1]
        suffixe = kmer[1::]
        graphe.add_edge(prefixe, suffixe, weight= kmer_dict[kmer])
    
    return graphe


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list :
        if not path :
            continue
        if delete_entry_node and delete_sink_node :
            graph.remove_nodes_from(path)
        elif delete_entry_node :
            graph.remove_nodes_from(path[:-1])
        elif delete_sink_node :
            graph.remove_nodes_from(path[1:])
        else :
            graph.remove_nodes_from(path[1:-1])
    return graph

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    if len(path_list) == 0:
        return graph

    meilleur_chemin = None
    ecart_type = statistics.stdev(weight_avg_list)

    if ecart_type > 0:
        index_meilleur_chemin = weight_avg_list.index(max(weight_avg_list))
        meilleur_chemin = path_list[index_meilleur_chemin]
    elif ecart_type == 0:
        ecart_type_longueur = statistics.stdev(path_length)

        if ecart_type_longueur > 0:
            index_meilleur_chemin = path_length.index(max(path_length))
            meilleur_chemin = path_list[index_meilleur_chemin]
        else:
            indice_aleatoire = random.randint(0, len(path_list) - 1)
            meilleur_chemin = path_list[indice_aleatoire]

    paths_to_remove = []
    for path in path_list:
        if path != meilleur_chemin:
            remove_path = path.copy()
            if not delete_entry_node:
                remove_path.pop(0)
            if not delete_sink_node:
                remove_path.pop(-1)
            paths_to_remove.extend(remove_path)
    
    graph.remove_nodes_from(paths_to_remove)

    return graph



def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    all_paths = list(nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node))
    path_lengths = [len(path) for path in all_paths]

    weight_avg_list = [
        sum(graph[path[i]][path[i + 1]].get('weight', 1) for i in range(len(path) - 1)) / max(len(path) - 1, 1)
        for path in all_paths
    ]

    return select_best_path(graph, all_paths, path_lengths, weight_avg_list)


def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble = False

    for node in graph.nodes:
        predecessors = list(graph.predecessors(node))

        if len(predecessors) > 1:
            for i in range(len(predecessors)):
                for j in range(i + 1, len(predecessors)):
                    ancestor_node = nx.lowest_common_ancestor(graph, predecessors[i], predecessors[j])
                    if ancestor_node is not None:
                        bubble = True
                        break
                if bubble:
                    break
        if bubble:
            break

    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, ancestor_node, node))

    return graph

def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    nodes_predecessors = [node for node, degree in graph.in_degree() if degree > 1]

    for node in nodes_predecessors:
        predecessors = list(graph.predecessors(node))

        if any(predecessor in starting_nodes for predecessor in predecessors):
            paths_start = [list(nx.all_simple_paths(graph, source=start, target=node))[0] for start in starting_nodes if nx.has_path(graph, start, node)]
            path_weights = [sum(graph[path[i]][path[i+1]].get('weight', 1) for i in range(len(path)-1)) for path in paths_start]
            max_weight = max(path_weights)
            main_paths = [path for i, path in enumerate(paths_start) if path_weights[i] == max_weight]

            for path in paths_start:
                if path not in main_paths:
                    for j in range(len(path) - 1):
                        graph.remove_edge(path[j], path[j+1])

    return graph

def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    for node in graph.nodes():
        successors = list(graph.successors(node))
        
        if len(successors) > 1:
            min_weight = min(graph[node][successor]['weight'] for successor in successors)
            graph.remove_edge(node, next(successor for successor in successors if graph[node][successor]['weight'] == min_weight))

    return graph

def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    list_nodes = []
    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            list_nodes.append(node)
    return(list_nodes)

def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    list_nodes = []
    for node in graph.nodes():
        if len(list(graph.successors(node))) == 0:
            list_nodes.append(node)
    return(list_nodes)

def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contig_len = []
    for starts_node in starting_nodes:
        for ending_node in ending_nodes:
            chemins = nx.all_simple_paths(graph, starts_node, ending_node)
            for chemin in chemins: # ["AT","TC","CT","TT"] -> ATCTT
                contig = chemin[0]
                for sequence in chemin[1:]:
                    contig += sequence[1]
                contig_len.append((contig, len(contig)))
    return contig_len

    
def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    with open(output_file, "w") as file :
        for i in range(len(contigs_list)):
            file.write(f">contig_{i} len={contigs_list[i][1]}\n")
            file.write(f"{textwrap.fill(contigs_list[i][0], width=80)}\n")


def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    file = args.fastq_file
    kmer = args.kmer_size

    print("Lecture du fichier et construction du graphe...")
    kmer_dict = build_kmer_dict(file, kmer)
    graph = build_graph(kmer_dict)
    #print("kmer_dict:", kmer_dict)

    print("Résolution des bulles...")
    graph = simplify_bubbles(graph)

    print("Résolution des pointes d'entrée et de sortie...")
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)

    print("Starting nodes:", starting_nodes)
    print("Ending nodes:", ending_nodes)

    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)

    print("Ecriture du/des contigs...")
    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs_list, args.output_file)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)


if __name__ == '__main__': # pragma: no cover
    main()
