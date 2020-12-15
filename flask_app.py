from flask import Flask, redirect
from flask import render_template
from flask import Markup
import pandas as pd
import json
from statistics import median

app = Flask(__name__)


@app.route('/')
def home():
    return redirect('/graph')


@app.route('/graph')
def graph():
    graph = open('templates/graph.html', 'r').read()
    return render_template('home.html', content=Markup(graph))


@app.route('/stats')
def stats():
    results = pd.read_csv('output/results.csv')
    results = results[results['read'] == 1.0]
    x = results['k'].to_list()
    y = results['n_cons'].to_list()
    ctgs = {'x': x, 'y': y, 'type': 'bar'}

    n50 = results['lens'].to_list()
    scores = []
    for a in n50:
        a = a.strip('[').strip(']')
        a = a.split(',')
        a = [int(n) for n in a]
        scores.append(median(a))

    als = results['al_scores'].to_list()
    asc = []
    for a in als:
        a = a.strip('[').strip(']').strip()
        a = a.split('.')
        l = []
        for n in a:
            c = n.strip('').strip('\n')
            if c != '':
                l.append(int(c))
        asc.append(sum(l) / len(l))

    n50 = {'x': x, 'y': scores}
    ascs = {'x': x, 'y': asc}

    stats = render_template('stats.html', contigs=ctgs, n50=n50, asc=ascs)
    return render_template('home.html', content=Markup(stats))


@app.route('/align')
def alignment():
    with open('output/k_30_th_2/best_alignment_1.json') as dt:
        data = json.load(dt)['align'][0]
    x, y = [], []
    for point in data:
        x.append(point[0])
        y.append(point[1])

    scatter = {
        'x': x,
        'y': y,
        'mode': 'markers',
        'type': 'scatter'
    }

    stats = render_template('alignment.html', scatter=scatter)
    return render_template('home.html', content=Markup(stats))


@app.route('/data')
def data():
    with open('output/k_30_th_2/path_file_0.json') as o:
        data1 = json.load(o)['path']
    with open('output_np/k_35_th_2/path_file_0.json') as o:
        data2 = json.load(o)['path']

    data = convert_graph_to_json([data1, data2])
    return data


def convert_graph_to_json(paths_list):
    nodes = set(paths_list[0] + paths_list[1])

    node_ref = {}
    for node in nodes:
        node_ref[node] = 0
        node_ref[node] = 0

    nodes = []
    for idx, key in enumerate(node_ref.keys()):
        nodes.append({
            "name": key,
            "label": key,
            "id": idx + 1
        })
        node_ref[key] = idx + 1

    edges_dic = {}
    edges = []
    for paths in paths_list:
        for i in range(len(paths) - 1):
            key = "s_{}_t_{}".format(node_ref[paths[i]], node_ref[paths[i + 1]])
            if key not in edges_dic.keys():
                edges_dic[key] = 1
            else:
                edges_dic[key] += 1

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
