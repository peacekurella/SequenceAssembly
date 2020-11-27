import argparse


def plot_de_bruijin_graph(graph):
    return None


def plot_discrepancies(mismatches):
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mismatches", default=None, help="mismatches file")
    parser.add_argument("--graph", default=None, help="de bruijin graph")