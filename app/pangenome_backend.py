import os
from os import walk
import re
import random

from abc import abstractmethod, ABCMeta
from collections import deque
from configparser import ConfigParser
from pathlib import Path
from typing import Optional, Union, Iterable

from fastapi import FastAPI, File, UploadFile
from fastapi.staticfiles import StaticFiles
from starlette.responses import HTMLResponse, FileResponse

import networkx
import pyfrost
from pyfrost import Kmer

GraphLike = Union[networkx.DiGraph, pyfrost.BifrostDiGraph]

app = FastAPI()
app.mount("/static/css", StaticFiles(directory="./client/build/static/css"), name="css")
app.mount("/static/js", StaticFiles(directory="./client/build/static/js"), name="js")

settings_file = os.environ.get("PGV_SETTINGS_FILE",
                               Path(__file__).with_name("settings.conf"))
app_config = ConfigParser()
app_config.read([settings_file])

data_path = 'graphs'


def kmerize_seq(k, seq: str, canonical: bool=False) -> Iterable[Kmer]:
    for pos in range(len(seq) - k + 1):
        kmer = Kmer(seq[pos:pos+k])

        if canonical:
            yield kmer.rep()
        else:
            yield kmer


class BaseGraphBackend(metaclass=ABCMeta):
    """Backend to obtain to load a pan-genome graph,
    and provide an API to obtain the graph neighborhood for
    a given query sequence.

    In the future we could provide backends for any kind of GFA graph,
    with an accompanying index."""

    @abstractmethod
    def get_neighborhood(self, query: str, radius: int=3) -> GraphLike:
        """
        Given a query sequence, return a NetworkX(-like) graph or subgraph object
        representing the graph neighborhood for that sequence.
        """


class PyfrostBackend(BaseGraphBackend):
    def __init__(self, graph_file):
        super().__init__()

        try:
            self.ccdbg = pyfrost.load(graph_file)

        except RuntimeError:
            self.ccdbg = None

    def get_graph(self):
        return self.ccdbg

    def get_neighborhood(self, query: str, radius: int=5) -> GraphLike:
        if not self.ccdbg:
            raise Exception("Could not load graph")

        neighborhood = networkx.DiGraph()

        for kmer in kmerize_seq(self.ccdbg.graph['k'], query):
            nodedict = self.ccdbg.find(kmer)
            if not nodedict:
                raise ValueError("K-mer does not exist in the graph")

            node = nodedict['head']

            if node in neighborhood:
                # Already added in a previous iteration
                continue

            # print(node)

            queue = deque()
            queue.append((0, node))

            while queue:
                level, node = queue.popleft()

                # print(level)
                # print(node)
                # print(type(node))
                # print(node.rep())
                # print(node.full_node())

                # ndata = self.ccdbg.nodes[node]

                neighborhood.add_node(node)

                if level < radius:
                    for succ in self.ccdbg.successors(node):
                        if succ in neighborhood:
                            continue

                        print("succ")
                        print(succ)

                        queue.append((level+1, succ))

                    for pred in self.ccdbg.predecessors(node):
                        if pred in neighborhood:
                            continue

                        print("pred")
                        print(pred)

                        queue.append((level+1, pred))

        # Add edges in the subgraph
        neighborhood.add_edges_from(self.ccdbg.out_edges(nbunch=list(neighborhood.nodes)))

        return neighborhood

    def get_nodes(self, data=True):
        if data:
            return (
                {"id": str(n), "sequence": data['unitig_sequence'], "tags": {"colors": list(data['colors'])}}
                for n, data in self.ccdbg.nodes(data=True)
            )
        else:
            return (
                {"id": str(n)} for n in self.ccdbg.nodes(data=False)
            )

    def get_links(self):
        return (
            {"source": str(u), "target": str(v)}
            for u, v in self.ccdbg.out_edges
        )

    def get_color_names(self):
        return self.ccdbg.graph['color_names']


def prepare_graphs():
    f = {}

    for (dirpath, dirnames, filenames) in walk(data_path):
        for file in filenames:
            if file.endswith(".gfa"):
                f[os.path.basename(dirpath)] = {'path': f'{dirpath}/{file}', 'graph': None}

    return f


def get_graph(graph_backends, graph_name):
    if graph_backends[graph_name]['graph'] is None:
        graph_backends[graph_name]['graph'] = PyfrostBackend(graph_backends[graph_name]['path'])

    return graph_backends[graph_name]['graph']


# Prepare to lazily load graph(s) in memory
g = prepare_graphs()


@app.get("/")
async def root():
    return FileResponse('client/build/index.html')


@app.get("/manifest.json")
async def manifest():
    return FileResponse('client/build/manifest.json')


@app.get("/logo192.png")
async def logo192():
    return FileResponse('client/build/logo192.png')


@app.get("/reference")
async def get_references():
    return {"references": list(app_config['references'].keys())}


@app.get("/reference/{color_id}")
async def get_reference(color_id: int):
    # TODO
    pass


@app.post("/graph/upload")
async def create_upload_file(file: UploadFile = File(...)):
    return {"filename": file.filename}


@app.get("/graph/test")
async def get_test_graph():
    return {
        "nodes": list(get_graph(g, 'small_test').get_nodes()),
        "links": list(get_graph(g, 'small_test').get_links()),
    }


@app.get("/graph/list")
async def list_graphs():
    return list(g.keys())


@app.get("/graph/describe/{graph_name}")
async def describe_graph(graph_name: str):
    if graph_name in g:
        gfa_file = g[graph_name]['path']
        bfg_colors = re.sub(".gfa", ".bfg_colors", gfa_file)
        is_cached = g[graph_name]['graph'] is not None

        return {'gfa': gfa_file, 'bfg': bfg_colors, 'is_cached': is_cached}

    return {'gfa': None, 'bfg': None, 'is_cached': None}


@app.get("/graph/nodes/{graph_name}")
async def get_nodes(graph_name: str):
    pb = get_graph(g, graph_name)

    return pb.get_nodes()


@app.get("/graph/links/{graph_name}")
async def get_links(graph_name: str):
    pb = get_graph(g, graph_name)

    return pb.get_links()


@app.get("/graph/colors/{graph_name}")
async def get_colors(graph_name: str):
    pb = get_graph(g, graph_name)

    return pb.get_color_names()


@app.get("/graph/random_node/{graph_name}")
async def get_random_node(graph_name: str):
    pb = get_graph(g, graph_name)

    nodes = list(pb.get_nodes(data=False))
    node = random.choice(nodes)

    return node


@app.get("/graph/neighborhood/{graph_name}")
async def get_neighborhood(graph_name: str, query: str, radius: int=5):
    pb = get_graph(g, graph_name)

    pb.get_neighborhood(query, radius)

    return 'hi'

