"""
pyComplexome - Proteome-scale inference of protein-protein interactions.

(C) Chris MacRaild, Monash University, 2022

A Bokeh app for exploring the complexome.
"""
from base64 import b64decode

import io
import os
import sys
from os.path import join, dirname

import numpy as np
import pandas as pd

import networkx as nx

from bokeh.models import ColumnDataSource, Circle, MultiLine, HoverTool, AutocompleteInput, DataTable, \
    TableColumn, Button, Column, Row, LabelSet, HTMLTemplateFormatter, CDSView, BooleanFilter, FileInput, CustomJS
from bokeh.plotting import figure, from_networkx
from bokeh.io import curdoc

from matplotlib.colors import LinearSegmentedColormap

from data import evidence_groups


score_palette = ['#9bb4e7', '#93a3e0', '#8b92d9', '#8381d1', '#7a71c9', '#7260c0', '#6a4fb6', '#623eac', '#5a2ca0']


def get_color_from_score(score):
    score_cm = LinearSegmentedColormap.from_list('score_cm', score_palette)
    return "#"+"".join([hex(int(256*c))[2:] for c in score_cm((score-0.6)/0.4)[:3]])


def get_score_formatter(my_col):
    template = f"""
        <div style="background:<%= 
            (function colorfromint(){{
                if ({my_col} < 0.6){{
                    return('{score_palette[0]}')}}
                else if ({my_col} < 0.65)
                    {{return('{score_palette[1]}')}}
                else if ({my_col} < 0.7)
                    {{return('{score_palette[2]}')}}
                else if ({my_col} < 0.75)
                    {{return('{score_palette[3]}')}}
                else if ({my_col} < 0.8)
                    {{return('{score_palette[4]}')}}
                else if ({my_col} < 0.85)
                    {{return('{score_palette[5]}')}}
                else if ({my_col} < 0.9)
                    {{return('{score_palette[6]}')}}
                else if ({my_col} < 0.95)
                    {{return('{score_palette[7]}')}}
                else
                    {{return('{score_palette[8]}')}}
                }}()) %>; 
            color: white;
            text-align: center"
            data-toggle="tooltip" 
            title="Overall confidence score for this interaction"> 
        <%= (value).toFixed(2) %>
        </div>
    """

    return HTMLTemplateFormatter(template=template)


def get_xlink_formatter(my_col):
    template = """
        <div style="background:<%= 
            (function colorfromint(){
                if (result_col == ''){
                    return('#ffffff')}
                else {
                    return('#37c837bb')}
                }()) %>; 
            color: white;
            text-align: center"> 
        <%= value %>
        </div>
    """.replace('result_col', my_col)

    return HTMLTemplateFormatter(template=template)


def get_evidence_formatter(my_col, hover_text):
    template = """
        <div style="background:<%= 
            (function colorfromint(){
                if (result_col < 0.1){
                    return('#fcea8cbb')}
                else if (result_col < 0.2)
                    {return('#ffd175bb')}
                else if (result_col < 0.3)
                    {return('#fdb860bb')}
                else if (result_col < 0.4)
                    {return('#f7a04cbb')}
                else if (result_col < 0.5)
                    {return('#ef893bbb')}
                else if (result_col < 0.6)
                    {return('#e4722cbb')}
                else if (result_col < 0.7)
                    {return('#d75b1ebb')}
                else if (result_col < 0.8)
                    {return('#c94312bb')}
                else if (result_col < 0.9)
                    {return('#ba2a06bb')}
                else if (result_col < 1.0)
                    {return('#a90000bb')}
                else
                    {return('#ffffff')}
                }()) %>; 
            color: white;
            text-align: center"
            data-toggle="tooltip" 
            title="hover"> 
        <%= (+value || 0).toFixed(2) %>
        </div>
    """.replace('result_col', my_col).replace('hover', hover_text)

    return HTMLTemplateFormatter(template=template)


evidence_descr = {'E': "Evers et al BNPage PCP dataset (https://doi.org/10.1038/s41467-021-23919-x)",
                    'H': "Hillier et al BNPage PCP dataset (https://doi.org/10.1016/j.celrep.2019.07.019)",
                    'M': 'MacRaild et al SEC PCP dataset (unpublished)',
                    'C': 'coexpression datasets'}
def get_hover_text(evidence_class):
    if evidence_class in evidence_descr:
        return f"""Percentile rank for this interaction in the {evidence_descr[evidence_class]}"""
    else:
        return f"""Percentile rank for this interaction in dataset {evidence_class}"""


def build_data_table(edge_ds, view=None, table=None):
    columns = []
    for col_name in ['protein1', 'protein2']:
        columns.append(TableColumn(field=col_name, title=col_name, width=90))
    for col_name in ['protein1_description', 'protein2_description']:
        columns.append(TableColumn(field=col_name, title=col_name))
    for col_name in edge_ds.column_names:
        if col_name == 'Crosslink':
            columns.append(TableColumn(field=col_name, title='Xlink', width=35,
                                       formatter=get_xlink_formatter(col_name)))
        elif col_name in evidence_groups:
            hover_text = get_hover_text(col_name)
            columns.append(TableColumn(field=col_name, title=col_name, width=35,
                                       formatter=get_evidence_formatter(col_name, hover_text)))
    columns.append(TableColumn(field='gbd_full_prob', title='Score', width=35,
                               formatter=get_score_formatter('gbd_full_prob')))

    # Calculate the available width for protein descriptions
    descr_width = (800 - (35*len(edge_ds.column_names) + 40))//2
    if descr_width < 150:
        descr_width = 150
    columns[2].width = descr_width
    columns[3].width = descr_width

    if table is None:
        if view is None:
            table = DataTable(source=edge_ds, columns=columns, width=800, index_position=None)
        else:
            table = DataTable(source=edge_ds, view=view, columns=columns, width=800, index_position=None)
    else:
        table.source = edge_ds
        table.columns = columns
        if view is not None:
            table.view = view

    table.autosize_mode = 'none'
    column_widths = sum([c.width for c in columns])
    if column_widths > 800:
        # If we need a horizontal scroll bar, keep the protein IDs fixed
        table.frozen_columns = 2
    #table.sizing_mode = 'fixed'

    table.background = 'white'

    return table


def filter_on_protein(attr, old, new):
    protein_mask = (df.protein1==new) | (df.protein2==new)
    update_table_on_selection(protein_mask)


def find_neighbours(event):
    idx = cds.selected.indices
    selected_proteins = set(cds.data['protein1'][idx]).union(cds.data['protein2'][idx])
    neighbour_mask = (df.protein1.isin(selected_proteins)) | (df.protein2.isin(selected_proteins))
    update_table_on_selection(neighbour_mask)


def update_table_on_selection(selected_mask):
    filtered_df = df.loc[selected_mask]
    neighbour_nodes = set(filtered_df.protein1).union(filtered_df.protein2)
    neighbour_edge_mask = df.protein1.isin(neighbour_nodes) & df.protein2.isin(neighbour_nodes)
    mask = selected_mask | neighbour_edge_mask
    view.filters = [BooleanFilter(mask), ]
    # ... and a hack workaround to force update [https://github.com/bokeh/bokeh/issues/11955]:
    table.visible = False
    table.visible = True


def remove_row(event):
    idx = cds.selected.indices
    mask = view.filters[0].booleans
    for i in idx:
        if not mask[i]: # debug
            print(f"Trying to remove missing row: index={idx}")
        mask[i] = False
    view.filters = [BooleanFilter(mask), ]
    # ... and a hack workaround to force update [https://github.com/bokeh/bokeh/issues/11955]:
    table.visible = False
    table.visible = True


def build_graph(cds, view):
    cp_graph = nx.Graph()

    selection_mask = view.filters[0].booleans
    indices = df.index[selection_mask]
    n_active_nodes = sum(selection_mask)
    # scale linewidth and node size according to the number of nodes to plot
    scale = 1 + 14/np.sqrt(n_active_nodes)

    node_descr = {}
    for idx in indices:
        n1 = cds.data['protein1'][idx]
        n2 = cds.data['protein2'][idx]
        node_descr[n1] = cds.data['protein1_description'][idx]
        node_descr[n2] = cds.data['protein2_description'][idx]
        edge_score = cds.data['gbd_full_prob'][idx]
        if ('Crosslink' in cds.data) and cds.data['Crosslink'][idx]:
            edge_xl = True
            edge_color = '#37c837'
            alpha = 0.6
        else:
            edge_xl = False
            edge_color = get_color_from_score(edge_score)
            alpha = 0.6
        line_width = scale
        cp_graph.add_edge(n1, n2, weight=((cds.data['gbd_full_prob'][idx] - 0.5)*2)**2,
                          edge_color=edge_color, alpha=alpha, line_width=line_width,
                          edge_score=edge_score, edge_xl=edge_xl)

    node_color = {}
    node_size = {}
    node_alpha = {}

    for node in cp_graph:
        node_alpha[node] = 1.0
        node_color[node] = 'lightgrey'
        node_size[node] = 5 * scale

    nx.set_node_attributes(cp_graph, node_descr, 'node_descr')
    nx.set_node_attributes(cp_graph, node_color, 'node_color')
    nx.set_node_attributes(cp_graph, node_size, 'node_size')
    nx.set_node_attributes(cp_graph, node_alpha, 'alpha')

    return cp_graph


def update_node_labels(x, y, text):
    """Replace existing node labels (if any) with labels text at position x,y

    x, y and text may be iterables of equal length"""

    label_source = ColumnDataSource(data={'x': x,
                                          'y': y,
                                          'label': text})
    labels = LabelSet(x='x', y='y', text='label', source=label_source)

    if len(plot.renderers) > 1:
        plot.renderers[1] = labels
    else:
        plot.renderers.append(labels)


def plot_graph(event):
    cp_graph = build_graph(cds, view)
    cp_graph_renderer = from_networkx(cp_graph, nx.spring_layout)
    cp_graph_renderer.node_renderer.glyph = Circle(size='node_size', fill_color='node_color', line_color='node_color',
                                                   line_alpha='alpha', fill_alpha='alpha')
    cp_graph_renderer.edge_renderer.glyph = MultiLine(line_color="edge_color", line_alpha='alpha',
                                                      line_width='line_width')

    plot.renderers = [cp_graph_renderer, ]

    x,y = zip(*cp_graph_renderer.layout_provider.graph_layout.values())
    n_nodes = len(cp_graph)
    labels = None

    if n_nodes < 15:
        # Graph is simple enough to label all nodes:
        update_node_labels(x, y, cp_graph_renderer.node_renderer.data_source.data['index'])

    elif protein_selector.value:
        # Label the node corresponding to the current value of the protein selector:
        try:
            idx = cp_graph_renderer.node_renderer.data_source.data['index'].index(protein_selector.value)
        except ValueError:
            # Selected protein not in plotted graph. Is this possible?
            pass
        update_node_labels([x[idx],], [y[idx],], [protein_selector.value,])

    elif cds.selected.indices:
        # Label nodes involved in selected edges:
        edge_idx = cds.selected.indices
        if len(edge_idx) <= 15:
            # Otherwise too many labels
            highlight_graph_edge(None, None, None)
            selected_proteins = set(cds.data['protein1'][edge_idx]).union(cds.data['protein2'][edge_idx])
            label_x = []
            label_y = []
            label_text = []
            for prot in selected_proteins:
                try:
                    idx = cp_graph_renderer.node_renderer.data_source.data['index'].index(prot)
                except ValueError:
                    pass
                else:
                    label_x.append(x[idx])
                    label_y.append(y[idx])
                    label_text.append(prot)
            update_node_labels(label_x, label_y, label_text)


def reset_table(event):
    view.filters[0] = BooleanFilter([True, ]*len(df))
    # ... and a hack workaround to force update [https://github.com/bokeh/bokeh/issues/11955]:
    table.visible = False
    table.visible = True
    #cds.data = dict(ColumnDataSource(df).data)


def highlight_graph_edge(attr, old, new):
    if plot.renderers:
        cp_graph_renderer = plot.renderers[0]
    else:
        # No plotted graph, so nothing to do
        return

    edge_data = dict(cp_graph_renderer.edge_renderer.data_source.data)
    node_data = dict(cp_graph_renderer.node_renderer.data_source.data)

    n_selected_edges = len(cds.selected.indices)
    if n_selected_edges < 15:
        # sufficently few selected edges, so label nodes:
        x,y = zip(*cp_graph_renderer.layout_provider.graph_layout.values())
        label_x = {}
        label_y = {}
        label_text = set()

    for cds_idx in cds.selected.indices:
        proteins = {cds.data['protein1'][cds_idx], cds.data['protein2'][cds_idx]}
        if n_selected_edges < 15:
            for prot in proteins:
                try:
                    prot_idx = node_data['index'].index(prot)
                except ValueError:
                    # Selected node is not currently plotted
                    continue
                label_x[prot] = x[prot_idx]
                label_y[prot] = y[prot_idx]
                label_text.add(prot)

        for edge_idx in range(len(edge_data['start'])):
            edge_proteins = {edge_data['start'][edge_idx], edge_data['end'][edge_idx]}
            if edge_proteins == proteins:
                edge_data['edge_color'][edge_idx] = 'orangered'
            elif ('Crosslink' in cds.data) and cds.data['Crosslink'][edge_idx]:
                edge_data['edge_color'][edge_idx] = '#37c837'
            else:
                edge_data['edge_color'][edge_idx] = get_color_from_score(edge_data['edge_score'][edge_idx])

    if n_selected_edges < 15:
        update_node_labels([label_x[k] for k in label_text],
                           [label_y[k] for k in label_text],
                           list(label_text))
    cp_graph_renderer.edge_renderer.data_source.data['edge_color'] = edge_data['edge_color']


protein_selector = AutocompleteInput(title="Filter on protein:", completions=[],
                                     min_characters=6, max_length=13, case_sensitive=False, )
protein_selector.on_change('value', filter_on_protein)

expand_neighbours = Button(label='Find neighbours', width=100, align=('start', 'end'))
expand_neighbours.on_click(find_neighbours)

remove_row_button = Button(label='Remove rows', width=100, align=('start', 'end'))
remove_row_button.on_click(remove_row)

reset_button = Button(label='Reset table', width=100, align=('start', 'end'))
reset_button.on_click(reset_table)

plot_graph_button = Button(label='Plot graph', width=100, align=('start', 'end'))
plot_graph_button.on_click(plot_graph)

node_hover_tool = HoverTool(tooltips=[("ID", "@index"), ("Description", "@node_descr")])

plot = figure()
#plot.add_tools(node_hover_tool, TapTool(), BoxSelectTool())
plot.add_tools(node_hover_tool)
plot.plot_width = 600
plot.plot_height = 600
plot.grid.grid_line_color = None
plot.xaxis.visible = False
plot.yaxis.visible = False


df = pd.DataFrame()
cds = ColumnDataSource(df)
view = CDSView(source=cds, filters=[BooleanFilter([True, ]*len(df))])
table = build_data_table(cds, view=view)


def load_data(attr, old, new):
    global df, cds, view, table
    if attr is None:
        # We have a filename:
        df = pd.read_csv(new, index_col=0, )
    else:
        # We have base64-encoded file contents from the FileInput widget
        decoded = b64decode(new)
        f = io.BytesIO(decoded)
        df = pd.read_csv(f, index_col=0, )
    df = df.reset_index()
    df.Crosslink = df.Crosslink.fillna('')
    cds = ColumnDataSource(df)
    cds.selected.on_change('indices', highlight_graph_edge)
    view = CDSView(source=cds, filters=[BooleanFilter([True, ] * len(df))])
    table = build_data_table(cds, view=view, table=table)
    proteins = sorted(set(df.protein1).union(df.protein2))
    protein_selector.completions = proteins


download_button = Button(label="Download .csv", button_type="success")
download_button.js_on_click(CustomJS(args=dict(table=table),
                                     code=open(join(dirname(__file__), "download.js")).read()))
controls = Row(protein_selector, expand_neighbours, remove_row_button, reset_button, plot_graph_button,
               sizing_mode='fixed', height=60, width=800)


# Load raw data:
if len(sys.argv) > 1:
    file_name = sys.argv[1]
else:
    file_name = join(dirname(__file__), 'data', 'gbt_draft_060622_Dnorm_edge_summary.csv')

if not os.path.exists(file_name):
    get_file = FileInput(accept='.csv')
    get_file.on_change('value', load_data)
    controls = Column(get_file, controls)
else:
    load_data(None, None, file_name)


left = Column(controls, table, download_button, width=800, height=600, sizing_mode='fixed')
curdoc().theme = 'light_minimal'
curdoc().add_root(Row(left, plot))