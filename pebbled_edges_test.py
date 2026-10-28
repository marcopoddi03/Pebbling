import re
import sys

def analyze_pebbled_edges(filename):
    edge_colors = {}

    with open(filename, "r") as file:
        for line in file:
            if not line.startswith(">"):
                continue

            fields = line[1:].split('|')

            edge_id = fields[0].strip()
            color = fields[2].strip()

            node_ids = re.search(r'>(.+?)([+-])>(.+?)([+-])', edge_id)
            if node_ids is None:
                continue

            edge = (
                node_ids.group(1),
                node_ids.group(2),
                node_ids.group(3),
                node_ids.group(4)
            )
            edge_colors.setdefault(edge, set()).add(color)

    edges_count = len(edge_colors)

    average_colors_per_edge = (
        sum(len(colors) for colors in edge_colors.values()) / edges_count
        if edges_count > 0
        else 0.0
    )

    return edges_count, average_colors_per_edge

filename = sys.argv[1]
edges_count, average_colors_per_edge = analyze_pebbled_edges(filename)

print(f"Number of pebbled edges: {edges_count}")
print(f"Average colors per edge: {average_colors_per_edge:.4f}")
