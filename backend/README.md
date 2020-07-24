# Pangenome visualization backend

Provides several API backends to obtain data for graph visualization and read alignments

We use [FastAPI](https://fastapi.tiangolo.com/) as API framework.

# Set up environment

```bash
python3 -mvenv venv
. venv/bin/activate
pip install numpy
pip install -r requirements.txt
```

# How to run

```bash
uvicorn pangenome_backend:app --reload
```

