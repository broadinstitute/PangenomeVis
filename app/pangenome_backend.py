import os
from abc import abstractmethod, ABCMeta
from collections import deque
from configparser import ConfigParser
from pathlib import Path
from typing import Optional, Union, Iterable

from fastapi import FastAPI
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
    def __init__(self):
        super().__init__()

        try:
            self.ccdbg = pyfrost.load(app_config['pyfrost']['graph'],
                                      app_config['pyfrost']['colors'])
        except RuntimeError:
            self.ccdbg = None

    def get_neighborhood(self, query: str, radius: int=5) -> GraphLike:
        if not self.ccdbg:
            raise Exception("Could not load graph")

        neighborhood = networkx.DiGraph()

        for kmer in kmerize_seq(self.ccdbg.graph['k'], query):
            node = self.ccdbg.find(kmer)
            if not node:
                raise ValueError("K-mer does not exist in the graph")

            if node in neighborhood:
                # Already added in a previous iteration
                continue

            queue = deque()
            queue.append((0, node))

            while queue:
                level, node = queue.popleft()

                ndata = self.ccdbg.nodes[node]
                neighborhood.add_node(node, **ndata)

                if level < radius:
                    for succ in self.ccdbg.successors(node):
                        if succ in neighborhood:
                            continue

                        queue.append((level+1, succ))

                    for pred in self.ccdbg.predecessors(node):
                        if pred in neighborhood:
                            continue

                        queue.append((level+1, pred))

        # Add edges in the subgraph
        neighborhood.add_edges_from(self.ccdbg.out_edges(nbunch=list(neighborhood.nodes)))

        return neighborhood

    def get_nodes(self):
        return (
            {"id": str(n), "sequence": data['unitig_sequence'], "tags": {"colors": list(data['colors'])}}
            for n, data in self.ccdbg.nodes(data=True)
        )

    def get_links(self):
        return (
            {"source": str(u), "target": str(v)}
            for u, v in self.ccdbg.out_edges
        )

    def get_color_names(self):
        return self.ccdbg.graph['color_names']


# Keep graph in memory
graph_backend = PyfrostBackend()


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


@app.get("/graph")
async def get_test_graph():
    return {
        "nodes": list(graph_backend.get_nodes()),
        "links": list(graph_backend.get_links()),
        "color_names": list(graph_backend.get_color_names())
    }