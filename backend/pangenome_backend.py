from configparser import ConfigParser
from pathlib import Path
from typing import Optional

from fastapi import FastAPI

app = FastAPI()


class Graph:
    def __init__(self):
        config_file = Path(__file__).with_name("graph.conf")
        self.config = ConfigParser()
        self.config.read([config_file])

    def build_graph():
        # Build compacted colored De Bruijn graph in memory
        pass


@app.get("/reference/{color_id}")
async def get_reference(color_id: Optional[int]):
    pass
