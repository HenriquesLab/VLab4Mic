import matplotlib.pyplot as plt
import seaborn as sns
from ..utils.transform.datatype import truncate


def pivot_dataframes(dataframe, axes_param_names):
    # extract individual dataframe per condition
    conditions = list(dataframe["Condition"].unique())
    subset_dataframes = [
        dataframe[dataframe["Condition"] == cond] for cond in conditions
    ]
    df_pivots = dict()
    # get mean and std accross parameter combinations of axes_param_names
    for condition_df in subset_dataframes:
        print(condition_df)
        condition_mean = condition_sd = None
        condition_mean = (
            condition_df.groupby(axes_param_names)["Metric"].mean().reset_index()
        )
        condition_sd = (
            condition_df.groupby(axes_param_names)["Metric"].std().reset_index()
        )
        # pivot dfs
        condition_mean_pivot = condition_mean.pivot(
            index="Fractional_defect", columns="Labelling_efficiency", values="Metric"
        ).round(4)
        condition_sd_pivot = condition_sd.pivot(
            index="Fractional_defect", columns="Labelling_efficiency", values="Metric"
        ).round(4)
        condition_name = condition_df["Condition"].unique()[0]
        print(condition_name)
        df_pivots[condition_name] = [condition_mean_pivot, condition_sd_pivot]
    return df_pivots


def sns_heatmap_pivots(
    df_pivots, titles = None, conditions_cmaps=None, annotations=False, cmaps_range="same", **kwargs
):
    conditions = list(df_pivots.keys())
    nconditions = len(conditions)
    annot_kws = {"size": 10, "rotation": 45}
    f, axes = plt.subplots(nconditions, 2, figsize=(12, 10))
    plot_num = 0
    if cmaps_range == "same":
        # min and max here correspond to SSIM
        hist_params = dict(vmin=0, vmax=1)
    elif cmaps_range == "each":
        hist_params = dict()
    if conditions_cmaps is None:
        conditions_cmaps = ["mako"] * n
    for n, cond in enumerate(conditions):
        print(cond, n)
        # mean
        sns.heatmap(
            df_pivots[cond][0],
            annot=annotations,
            annot_kws=annot_kws,
            ax=axes[n, 0],
            cmap=conditions_cmaps[n],
            xticklabels=df_pivots[cond][0].columns.values.round(3),
            yticklabels=df_pivots[cond][0].index.values.round(3),
            **hist_params,
        )
        axes[n, 0].set_title(titles["category"]+ ": " + cond + ". Mean Metric")
        # std
        sns.heatmap(
            df_pivots[cond][1],
            annot=annotations,
            annot_kws=annot_kws,
            ax=axes[n, 1],
            cmap=conditions_cmaps[n],
            xticklabels=df_pivots[cond][1].columns.values.round(3),
            yticklabels=df_pivots[cond][1].index.values.round(3),
        )
        axes[n, 1].set_title(titles["category"]+ ": " + cond + ". Std Dev Metric")
    f.tight_layout()


def show_references(references):
    n_conditions = len(list(references.keys()))
    f, axes = plt.subplots(1, n_conditions, figsize=(12, 10))
    i = 0
    for cond, img in references.items():
        axes[i].imshow(img, cmap="grey")
        axes[i].set_title(f"Reference for: {cond}")
        i = i + 1


def show_example_test(queries, params, condition="STED_demo", replica_number=1, query_variant=0):    #
    param_values = [truncate(p, 6) for p in params]
    print(param_values)
    combination_pars = [str(val) for val in param_values]
    print(combination_pars)
    combination_name = ",".join(combination_pars)
    print(combination_name)
    plt.imshow(queries[combination_name][replica_number][condition][query_variant], cmap="grey")
    plt.title(combination_name)
