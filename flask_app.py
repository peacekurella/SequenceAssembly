from flask import Flask
from flask import render_template
from flask import Markup
import json

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('home.html', content=None)


@app.route('/graph')
def graph():
    graph = open('templates/graph.html', 'r').read()
    return render_template('home.html', content=Markup(graph))


@app.route('/stats')
def stats():
    stats = open('templates/stats.html', 'r').read()
    return render_template('home.html', content=Markup(stats))


@app.route('/align')
def alignment():
    align = open('templates/graph.html', 'r').read()
    return render_template('home.html', content=Markup(align))


@app.route('/data')
def data():
    with open('output/k_30_th_2/graph.json') as o:
        data = json.load(o)
    return data


def convert_graph_to_json(edges):
    node_ref = {}
    for edge in edges:
        node_ref[edge[:-1]] = 0
        node_ref[edge[1:]] = 0

    nodes = []
    for idx, key in enumerate(node_ref.keys()):
        nodes.append({
            "name": key,
            "label": key,
            "id": idx + 1
        })
        node_ref[key] = idx + 1

    edges_dic = {}
    for edge in edges:
        key = "s_{}_t_{}".format(node_ref[edge[:-1]], node_ref[edge[1:]])
        if key not in edges_dic.keys():
            edges_dic[key] = 1
        else:
            edges_dic[key] += 1

    edges = []
    for key in edges_dic.keys():
        s, t = key.split("_")[1], key.split("_")[3],
        edges.append({
            "source": s,
            "target": t,
            "type": edges_dic[key]
        })

    final_op = {
        "nodes": nodes,
        "links": edges
    }

    return final_op


if __name__ == "__main__":
    app.run(debug=True)
