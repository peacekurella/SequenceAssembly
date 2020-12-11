from flask import Flask
from flask import render_template
from flask import Markup

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('home.html', user_name='Hi!', content=None)


@app.route('/stats')
def stats():
    stats = open('templates/graph.html', 'r').read()
    return render_template('home.html', user_name='Hi!', content=Markup(stats))


@app.route('/data')
def data():
    data = {
        "nodes": [
            {
                "name": "Peter",
                "label": "Person",
                "id": 1
            },
            {
                "name": "Michael",
                "label": "Person",
                "id": 2
            },
            {
                "name": "Neo4j",
                "label": "Database",
                "id": 3
            },
            {
                "name": "Graph Database",
                "label": "Database",
                "id": 4
            }
        ],
        "links": [
            {
                "source": 1,
                "target": 2,
                "type": "KNOWS",
                "since": 2010
            },
            {
                "source": 1,
                "target": 3,
                "type": "FOUNDED"
            },
            {
                "source": 2,
                "target": 3,
                "type": "WORKS_ON"
            },
            {
                "source": 3,
                "target": 4,
                "type": "IS_A"
            }
        ]
    }
    return data


if __name__ == "__main__":
    app.run(debug=True)