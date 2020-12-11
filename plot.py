import dash
import dash_core_components as dcc
import dash_html_components as html
import networkx as nx
import plotly.graph_objs as go


if __name__ == "__main__":
    edges = [
        (1, 2),
        (2, 3),
        (3, 2),
        (3, 4),
        (4, 2),
        (2, 3)
    ]

    nodes = [
        1,
        2,
        3,
        4,
        5
    ]

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    node_positions = nx.kamada_kawai_layout(G)

    xn, yn = [], []
    for node in node_positions:
        x, y = node_positions[node]
        xn.append(x)
        yn.append(y)

    xe, ye = [], []
    for edge in G.edges:
        x0, y0 = node_positions[edge[0]]
        x1, y1 = node_positions[edge[1]]
        xe.extend([x0, x1, None])
        ye.extend([y0, y1, None])

    edge_args = {
        'x': xe,
        'y': ye,
        'mode': 'lines',
        'line': {'width': 2},
        'hoverinfo': 'none',
        'name': 'lines'
    }
    edge_trace = go.Scatter(**edge_args)

    node_args = {
        'x': xn,
        'y': yn,
        'mode': 'markers',
        'name': 'k-mers',
        'marker': {'symbol': 'circle', 'size': 12},
        'text': nodes,
        'hoverinfo': 'text'
    }
    node_trace = go.Scatter(**node_args)

    axis = {
        'showline': False,
        'showgrid': False,
        'showticklabels': False,
    }

    layout = {
        'title': 'de-Bruijin graph',
        'width': 900,
        'height': 600,
        'showlegend': False,
        'margin': {'t': 100},
        'hovermode': 'closest',
        'plot_bgcolor': 'rgba(0,0,0,0)',
    }

    layout = go.Layout(**layout)

    data = [
        node_trace,
        edge_trace
    ]
    fig = go.Figure(data, layout)
    fig.update_xaxes(**axis)
    fig.update_yaxes(**axis)

    app = dash.Dash(__name__)
    app.layout = html.Div([
        html.Div([
            dcc.Graph(id="graph-1", figure=fig),
            dcc.Graph(id="graph-2", figure=fig),
        ])
    ])
    app.run_server(debug=True)
