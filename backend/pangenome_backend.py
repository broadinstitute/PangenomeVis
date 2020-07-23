from abc import abstractmethod, ABCMeta
from collections import deque
from configparser import ConfigParser
from pathlib import Path
from typing import Optional, Union, Iterable

import networkx
import pyfrost
from pyfrost import Kmer

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

GraphLike = Union[networkx.DiGraph, pyfrost.BifrostDiGraph]

app = FastAPI()
app.mount("/static", StaticFiles(directory="static"), name="static")

app_config = ConfigParser()
app_config.read([Path(__file__).with_name("settings.conf")])


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

        self.ccdbg = pyfrost.load(app_config['pyfrost']['graph'],
                                  app_config['pyfrost']['colors'])

    def get_neighborhood(self, query: str, radius: int=5) -> GraphLike:
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


# Keep graph in memory
graph_backed = PyfrostBackend()


@app.get("/reference/{color_id}")
async def get_reference(color_id: Optional[int], Depe):
    if color_id:
        # TODO
        pass
    else:
        return [app_config['references'].keys()]
