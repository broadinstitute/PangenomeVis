# Pangenome visualization backend

Provides several API backends to obtain data for graph visualization and read alignments

We use [FastAPI](https://fastapi.tiangolo.com/) as API framework.

# Set up environment

```bash
python3 -mvenv venv
. venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

# How to run locally

```bash
./deploy_locally.sh
```

# How to run on Google App Engine (must have Broad Institute credentials)

```bash
./deploy_gae.sh
```
