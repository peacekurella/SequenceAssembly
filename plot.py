import argparse
import igraph as ig
import plotly.graph_objs as go


def plot_de_bruijin_graph(graph):
    nodes = [
        {'name': '1', 'group': 1},
        {'name': '2', 'group': 2},
        {'name': '3', 'group': 3},
        {'name': '4', 'group': 1},
        {'name': '5', 'group': 2},
    ]

    edges = [
        (1, 5),
        (2, 5),
        (5, 1),
        (4, 3),
        (2, 1)
    ]

    group = [
        1,
        2,
        3
    ]

    labels = [
        1,
        2,
        3,
        4,
        5
    ]

    G = ig.Graph(edges, directed=True)
    layt = G.layout('kk', dim=3)

    xn = [layt[n][0] for n in range(len(layt))]
    yn = [layt[n][1] for n in range(len(layt))]
    zn = [layt[n][2] for n in range(len(layt))]

    xe = [[layt[e[0]][0], layt[e[1]][0], None] for e in edges]
    ye = [[layt[e[0]][1], layt[e[1]][1], None] for e in edges]
    ze = [[layt[e[0]][2], layt[e[1]][2], None] for e in edges]

    edge_args = {
        'x': xe,
        'y': ye,
        'z': ze,
        'mode': 'lines',
        'line': {'color': 'rgb(125, 125, 125)', 'width': 3},
        'hoverinfo': 'none',
        'name': 'lines'
    }
    edge_trace = go.Scatter3d(**edge_args)

    node_args = {
        'x': xn,
        'y': yn,
        'z': zn,
        'mode': 'markers',
        'name': 'k-mers',
        'marker': {'symbol': 'circle', 'size': 6, 'color': group, 'line': {'color': 'rgb(50, 50, 50)', 'width': 0.5}},
        'text': labels,
        'hoverinfo': 'text'
    }
    node_trace = go.Scatter3d(**node_args)

    layout = {
        'title': 'de-Bruijin graph',
        'width': 1000,
        'height': 1000,
        'showlegend': False,
        'margin': {'t': 100},
        'hovermode': 'closest'
    }

    layout = go.Layout(**layout)

    data = [
        node_trace,
        edge_trace
    ]
    fig = go.Figure(data, layout)

    fig.show()


def plot_discrepancies(mismatches):
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mismatches", default=None, help="mismatches file")
    parser.add_argument("--graph", default=None, help="de bruijin graph")

    plot_de_bruijin_graph(None)
