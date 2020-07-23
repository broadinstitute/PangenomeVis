from flask import Flask, render_template

app = Flask(__name__)


@app.route("/")
def home():
    return render_template("index.html")


@app.route("/review/private")
def review_private():
    return "Hello, Review (private)"


@app.route("/review/public")
def review_public():
    return render_template("review_public.html")


if __name__ == "__main__":
    app.run(debug=True)
