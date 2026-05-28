import yaml
import os
import sys
import importlib
import pandas as pd

# Ensure the parent directory is in the path so we can import from models
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from utils.plot_utils import Plotter

def load_config(filepath="config.yaml"):
    with open(filepath, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def run_stage(stage_name, plot_configs, input_data=None):
    """
    通用化的模組執行器。
    自動根據 stage_name 去尋找對應的資料夾與模組。
    可以接收前一個 stage 的輸出 (input_data) 作為輸入。
    """
    print(f"\n--- Running Stage: {stage_name} ---")
    config = load_config()
    stage_config = config.get(stage_name, {})
    
    if not stage_config:
        print(f"Warning: config.yaml 中找不到 {stage_name} 的設定，跳過執行。")
        return None
        
    # 動態載入該階段的 process.py 模組
    try:
        process_module = importlib.import_module(f"models.{stage_name}.process")
    except ImportError as e:
        print(f"Error: 找不到 {stage_name} 的 process.py 模組。({e})")
        return None

    # 執行運算，傳入設定檔與前一階段的數據
    result_df = process_module.run(stage_config, input_data=input_data)
    
    # 儲存中間數據
    intermediate_dir = os.path.join("data", "02_intermediate")
    os.makedirs(intermediate_dir, exist_ok=True)
    intermediate_path = os.path.join(intermediate_dir, f"{stage_name}_result.csv")
    result_df.to_csv(intermediate_path, index=False)
    print(f"Saved intermediate data to: {intermediate_path}")
    
    # 產生圖表
    print(f"\n--- Generating Plots for {stage_name} ---")
    plotter = Plotter(output_dir=os.path.join("outputs", f"{stage_name}_plots"))
    
    for plot_cfg in plot_configs:
        plot_type = plot_cfg.get("type")
        
        # If the result dataframe has an 'Analysis_Type' column, filter it first
        plot_data = result_df
        if 'analysis_filter' in plot_cfg and 'Analysis_Type' in result_df.columns:
            plot_data = result_df[result_df['Analysis_Type'] == plot_cfg['analysis_filter']]
            
        if plot_type == "line":
            plotter.plot_line(
                data=plot_data, 
                x_col=plot_cfg.get("x_col", "Time"),
                y_col=plot_cfg["y_col"],
                hue_col=plot_cfg["hue_col"],
                title=plot_cfg["title"], 
                filename=plot_cfg["filename"],
                legend_title=plot_cfg.get("legend_title"),
                x_label=plot_cfg.get("x_label"),
                log_x=plot_cfg.get("log_x", False)
            )
        elif plot_type == "dose_response":
            plotter.plot_dose_response(
                data=plot_data, 
                x_col=plot_cfg["x_col"],
                y_col=plot_cfg["y_col"], 
                title=plot_cfg["title"], 
                filename=plot_cfg["filename"],
                x_label=plot_cfg.get("x_label"),
                log_x=plot_cfg.get("log_x", False)
            )
            
    # 回傳結果給下一個 model 使用
    return result_df

def run_pipeline(models_to_run):
    print("--- Pipeline Execution Started ---")
    
    result_0_5 = None
    if "model_0_5" in models_to_run or "all" in models_to_run:
        # 執行 Model 0.5
        model_0_5_plots = [
            {"type": "line", "x_col": "Time", "y_col": "Cas12a", "hue_col": "Ara_Concentration", "title": "Cas12a Expression Dynamics over Time", "filename": "cas12a_time_series.png", "legend_title": "[Ara]", "x_label": "Time (minutes)"},
            {"type": "line", "x_col": "Time", "y_col": "mRNA", "hue_col": "Ara_Concentration", "title": "Cas12a mRNA Dynamics over Time", "filename": "mRNA_time_series.png", "legend_title": "[Ara]", "x_label": "Time (minutes)"},
            {"type": "dose_response", "x_col": "Ara_Concentration", "y_col": "Cas12a", "title": "Steady State Cas12a vs Arabinose Concentration", "filename": "cas12a_dose_response.png", "x_label": "Arabinose Concentration [Ara]"}
        ]
        result_0_5 = run_stage("model_0_5", model_0_5_plots)
    else:
        path = os.path.join("data", "02_intermediate", "model_0_5_result.csv")
        if os.path.exists(path):
            result_0_5 = pd.read_csv(path)
            
    result_1_0 = None
    if "model_1_0" in models_to_run or "all" in models_to_run:
        if result_0_5 is not None:
            model_1_0_plots = [
                {
                    "type": "line", 
                    "analysis_filter": "Activation",
                    "x_col": "Target_DNA", 
                    "y_col": "Activated_Cas12a", 
                    "hue_col": "Ara_Concentration", 
                    "title": "Cas12a Activation vs Target DNA", 
                    "filename": "activation_curve.png", 
                    "legend_title": "Inducer [Ara]",
                    "x_label": "Target DNA (nM)",
                    "log_x": True
                },
                {
                    "type": "line", 
                    "analysis_filter": "Cleavage",
                    "x_col": "ssDNA", 
                    "y_col": "Cleavage_Rate", 
                    "hue_col": "Target_DNA", 
                    "title": "Cleavage Rate vs ssDNA Substrate", 
                    "filename": "cleavage_kinetics.png", 
                    "legend_title": "Target DNA (nM)",
                    "x_label": "ssDNA Substrate (nM)",
                    "log_x": True
                },
                {
                    "type": "dose_response", 
                    "analysis_filter": "LOD",
                    "x_col": "Ara_Concentration", 
                    "y_col": "LOD", 
                    "title": "Limit of Detection (LOD) vs Inducer Concentration", 
                    "filename": "lod_curve.png", 
                    "x_label": "Arabinose Concentration [Ara]",
                    "log_x": False
                }
            ]
            result_1_0 = run_stage("model_1_0", model_1_0_plots, input_data=result_0_5)
    else:
        path = os.path.join("data", "02_intermediate", "model_1_0_result.csv")
        if os.path.exists(path):
            result_1_0 = pd.read_csv(path)

    if "model_2" in models_to_run or "all" in models_to_run:
        model_2_plots = [
            {
                "type": "line", 
                "x_col": "Time", 
                "y_col": "ssDNA", 
                "hue_col": "Target_DNA", 
                "title": "ssDNA Accumulation via RCR", 
                "filename": "rcr_ssDNA_time_series.png", 
                "legend_title": "Target DNA (nM)",
                "x_label": "Time (minutes)"
            }
        ]
        run_stage("model_2", model_2_plots, input_data=result_1_0)

    if "model_3" in models_to_run or "all" in models_to_run:
        model_3_plots = [
            {
                "type": "line", 
                "x_col": "Time", 
                "y_col": "Plasmid_Concentration", 
                "hue_col": "Target_DNA", 
                "title": "Sensor Plasmid Degradation over Time", 
                "filename": "plasmid_degradation_time_series.png", 
                "legend_title": "Target DNA (nM)",
                "x_label": "Time (minutes)"
            }
        ]
        run_stage("model_3", model_3_plots, input_data=result_1_0)
    
    print("\n--- Pipeline Execution Completed! ---")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="iGEM Pipeline Orchestrator")
    parser.add_argument("--models", nargs="+", default=["all"], help="Models to run, e.g., model_0_5 model_1_0 model_2 model_3 all")
    args = parser.parse_args()
    run_pipeline(args.models)
