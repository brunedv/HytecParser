import matplotlib.pyplot as plt  # plt.FUNCTIONNAE
import numpy as np  # outils numériques
import pandas as pd


def find_header(file_path):
    file1 = open(file_path, "r")

    current_line = file1.readline()
    row_number = 0
    while not ("sample Nr." in current_line):
        current_line = file1.readline()
        row_number += 1
    data_header = pd.read_csv(file_path, header=None, nrows=row_number - 1)

    list_col_name = []
    list_unit = []

    if ".flc" in file_path:
        list_col_name.append("Loc")
        list_unit.append(None)
    for j in range(len(data_header)):
        if data_header.loc[j, 0].find("column") > 0:
            col_name_unit = data_header.loc[j, 0].split(":")[1].split()
            if len(col_name_unit) < 2:
                list_col_name.append(col_name_unit[0])
                list_unit.append(None)
            else:
                list_col_name.append(col_name_unit[0])
                list_unit.append(col_name_unit[1])
    df_col_unit = pd.DataFrame({"column_name": list_col_name, "unit": list_unit})

    index_duplicates = df_col_unit.duplicated(subset=["column_name"])
    df_col_unit.loc[index_duplicates, "column_name"] = (
        df_col_unit.loc[index_duplicates, "column_name"]
        + df_col_unit.loc[index_duplicates, "unit"]
    )

    return df_col_unit, row_number


def read_data(file_path):
    df_col_unit, row_number = find_header(file_path)

    data_brut = pd.read_csv(
        file_path,
        sep="\s+",
        names=df_col_unit.column_name.values,
        skiprows=row_number - 1,
    )
    data_brut["time"] = data_brut[data_brut.iloc[:, 1] == "time"].iloc[:, 3]
    data_brut["time"] = data_brut["time"].fillna(method="ffill")
    data_brut.dropna(inplace=True)
    list_bad_col = data_brut.columns[data_brut.dtypes == "object"]
    if len(list_bad_col) > 1:
        data_brut[list_bad_col] = data_brut[list_bad_col].astype(float)
    return data_brut, df_col_unit, row_number


class plot_var:

    logscale: bool
    fac: float
    acidity: bool
    setxmin: bool
    setymin: bool
    setxmax: bool
    setymax: bool
    settmin: bool
    settmax: bool
    tmin: float
    tmax: float
    xmin: float
    ymin: float
    xmax: float
    ymax: float  # borner les axes xmin, xmax, ymin, ymax, tmin, tmax
    color: str
    name: str  # name: comme ça apparait dans hytec et les fichiers .res
    fancy_name: str  # pour renommer de manière "fancy"
    fname: str  # pour le nom du fichier qui va être généré
    transient: bool  # est-ce que tu veux que ce soit à un endroit au cour sdu temps
    profile: bool  # est-ce que veux profil spatial après un an
    to_plot: bool  # option inutile:
    unit: str  # l'unité que tu peux spécifier de la grandeur (ex: mg/l)
    type_sort: int

    def __init__(
        self,
        name,
        unit="",
        fancy_name=" ",
        fname=" ",
        to_plot=True,
        logscale=False,
        fac=1,
        setxmin=False,
        acidity=False,
        setymin=False,
        setxmax=False,
        setymax=False,
        transient=True,
        type_sort=0,
        profile=True,
        xmin=-1,
        ymin=-1,
        xmax=-1,
        ymax=-1,
        settmin=False,
        settmax=False,
        tmin=-1,
        tmax=-1,
        color="",
    ):
        self.logscale = logscale
        self.fac = fac
        self.setxmin = setxmin
        self.setymin = setymin
        self.setxmax = setxmax
        self.setymax = setymax
        self.settmin = settmin

        self.settmax = settmax
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax
        self.tmin = tmin
        self.tmax = tmax
        self.color = color
        self.type_sort = type_sort
        self.to_plot = to_plot
        if self.xmin != -1:

            self.setxmin = True
        if self.ymin != -1:
            self.setymin = True
        if self.xmax != -1:
            self.setxmax = True
        if self.ymax != -1:
            self.setymax = True
        if self.tmax != -1:
            self.settmax = True
        if self.tmin != -1:
            self.settmin = True
        self.acidity = acidity
        self.name = name
        self.unit = unit
        if fancy_name == " ":
            self.fancy_name = self.name
        else:
            self.fancy_name = fancy_name
        if fname == " ":
            self.fname = self.name
        else:
            self.fname = fname
        self.transient = transient
        self.profile = profile


class Simulation:
    fold: str
    legend_label: str
    plot_time: float
    color: str

    def __init__(self, fold, legend_label=" ", plot_time=[], color=""):
        self.fold = fold
        self.legend_label = legend_label
        self.plot_time = []
        self.color = color
        if legend_label == " ":
            self.legend_label = self.fold


class Combined_plot:
    title: str
    logscale: bool
    ylab: str
    fname: str
    setxmin: bool
    setymin: bool
    setxmax: bool
    setymax: bool
    xmin: float
    ymin: float
    xmax: float
    ymax: float

    def __init__(
        self,
        title,
        fname,
        ylab,
        logscale=False,
        setxmin=False,
        setymin=False,
        setxmax=False,
        setymax=False,
        xmin=-1,
        ymin=-1,
        xmax=-1,
        ymax=-1,
    ):
        self.title = title
        self.fname = fname
        self.ylab = ylab
        self.logscale = logscale
        self.setxmin = setxmin
        self.setymin = setymin
        self.setxmax = setxmax
        self.setymax = setymax
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax
        if self.xmin != -1:
            self.setxmin = True
        if self.ymin != -1:
            self.setymin = True
        if self.xmax != -1:
            self.setxmax = True
        if self.ymax != -1:
            self.setymax = True
