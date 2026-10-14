#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:06:46 2023
# Running python scripts in the linux background :
nohup python -u recycle_ISMC_wrfrst.py > recycle_ISMC_wrfrst.out 2>&1 &
tail -f recycle_ISMC_wrfrst.out
@author: veiga5
"""
from netCDF4 import Dataset
import os
import subprocess
import time
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import re
from scipy.interpolate import griddata, RegularGridInterpolator
from time import sleep

# 读取静态数据
def read_geo_info(DOMAIN_PATH):
    dst_nc = Dataset(os.path.join(DOMAIN_PATH, "geo_em.d01.nc"))
    dst_lat = dst_nc.variables['XLAT_M'][0,:,:]
    dst_lon = dst_nc.variables['XLONG_M'][0,:,:]
    landmask = dst_nc.variables['LANDMASK'][0,:,:]
    dst_nc.close()
    # 读取smcmax
    soil_nc = Dataset(os.path.join(DOMAIN_PATH, "soil_properties.nc"))
    smcmax = soil_nc.variables['smcmax'][0, 0, :, :]
    soil_nc.close()
    return dst_lat, dst_lon, landmask, smcmax 

def replace_zeros(array, missing_value):
    rows, cols = array.shape
    for i in range(rows):
        for j in range(cols):
            if array[i, j] == missing_value:
                # Find non-zero neighbors
                neighbors = []
                if i > 0 and array[i-1, j] != missing_value:
                    neighbors.append(array[i-1, j])
                if i < rows-1 and array[i+1, j] != missing_value:
                    neighbors.append(array[i+1, j])
                if j > 0 and array[i, j-1] != missing_value:
                    neighbors.append(array[i, j-1])
                if j < cols-1 and array[i, j+1] != missing_value:
                    neighbors.append(array[i, j+1])
                # Replace zero with a random non-zero neighbor
                if neighbors:
                    array[i, j] = np.random.choice(neighbors)
                else:
                    # If no non-zero neighbors found, replace with a random value
                    array[i, j] = np.random.randint(0.1, 0.5)  # Adjust range as needed
    return array

# 读取原始土壤湿度数据
def Read_Src_SM_data(SM_PATH, date, SM_data):
    lat_name = "lat"
    lon_name = "lon"
    if SM_data == 'CSMX':
        VAR_NAME = 'SM'
        sm_file = os.path.join(SM_PATH, f"{date[0:8]}.CSMX_000010.nc")

    with Dataset(sm_file) as src_nc:
        src_sm = src_nc.variables[f'{VAR_NAME}'][:]
        src_lat = src_nc.variables[f'{lat_name}'][:]
        src_lon = src_nc.variables[f'{lon_name}'][:]
    return src_sm, src_lat, src_lon

def interpolate_by_griddata(src_sm, src_lat, src_lon, dst_points, scale_factor, REGRIDmethod, 
           dst_lat, dst_lon, landmask, smcmax):
    print(f"Start interpolate by {REGRIDmethod} griddata method")
    StartInterp = time.time()
    src_lon_grid, src_lat_grid = np.meshgrid(src_lon, src_lat)
    src_points = np.array([src_lon_grid.flatten(), src_lat_grid.flatten()]).T
    # print(src_points.shape, src_values.shape, dst_points.shape)
    dst_sm_ori = griddata(src_points, (src_sm).flatten(), dst_points, method=REGRIDmethod)
    dst_sm_ori = dst_sm_ori.reshape((dst_lat.shape[0], dst_lon.shape[1]))
    dst_sm_ori *= scale_factor

    dst_sm_ori[landmask==0] = 1 # 海洋赋值为1
    print(f"Interpolate Data Elapse {time.time()-StartInterp} seconds")
    return dst_sm_ori
 
def write_SMdata_to_wrfrst(DST_SM_ORI, RUN_PATH, formatted_date, smcmax):
    # Copy file
    wrfrst_nc =  os.path.join(
        RUN_PATH,
        f"wrfrst_d01_{formatted_date}_00:00:00"
    )
    with Dataset(wrfrst_nc, 'r+') as wrfrstOutnc:
        print(wrfrstOutnc.variables['SMOIS'].shape, DEPTH_ID)
        varCtrl = wrfrstOutnc.variables['SMOIS'][0,DEPTH_ID,:,:] 
        DST_SM = DST_SM_ORI[DEPTH_ID,:,:]
        #---------------------------根据smcmax进行质量控制----------------------
        DST_SM = np.where(DST_SM > smcmax, smcmax, DST_SM)
        # ---------------------------替换负值为原始值-----------------------------
        DST_SM[DST_SM <= 0] = varCtrl[DST_SM <= 0]
        varCtrl = np.where(varCtrl > smcmax, smcmax, varCtrl)
        # ---------------------------替换初始土壤湿度-----------------------------
        wrfrstOutnc.variables['SMOIS'][0,DEPTH_ID,:,:] = DST_SM
        wrfrstOutnc.variables['SH2O'][0,DEPTH_ID,:,:] = DST_SM
    pass

def Replace_SM_data(RUN_PATH, SM_PATH, SM_data, date_str, dst_lon, dst_lat, landmask, smcmax):
    print("Start Replace SM data")
    # -------------------- 读取原始土壤湿度数据 -------------------------
    src_sm, src_lat, src_lon = Read_Src_SM_data(SM_PATH, date_str, SM_data)
    # ------------------- 插值到WRF经纬度网格 ------------------------
    DEPTH_ID =0 # 替换土壤湿度的垂直层次（1~4层）
    x_dim, y_dim = int(dst_lon.shape[1]), int(dst_lat.shape[0])
    z_dim =4
    scale_factor = 1
    REGRIDmethod = 'nearest'
    # 提前计算目标网格点阵列
    dst_points = np.array([dst_lon.flatten(), dst_lat.flatten()]).T
    # 强制类型转换
    dst_points = dst_points.astype(np.float32)
    DST_SM_ORI = np.zeros(( z_dim, y_dim, x_dim), dtype=np.float64)
    DST_SM_ORI[DEPTH_ID,:] = interpolate_by_griddata(src_sm, src_lat, src_lon, dst_points, scale_factor, REGRIDmethod, dst_lat, dst_lon, landmask, smcmax)
    # 使用 datetime 解析字符串并直接格式化输出
    formatted_date = datetime.strptime(date_str[0:8], "%Y%m%d").strftime("%Y-%m-%d")
    # --------------------- 将插值好的土壤湿度数据写入wrfrst文件 -----------------------
    print("Start write SM data to wrfrst file ...")
    write_SMdata_to_wrfrst(DST_SM_ORI, RUN_PATH, formatted_date, smcmax)
    pass

def replace_namelist(RUN_PATH, date_str, Replace_Interval_Hours):
    # 文件路径
    namelist_input_file_path = os.path.join(RUN_PATH, 'namelist.input')

    # 计算开始时间 
    start_year  = date_str[0:4]
    start_month = date_str[4:6]
    start_day   = date_str[6:8]	
    start_hour = date_str[8:10]	

    # 计算结束时间 
    start_dt = datetime(int(start_year), int(start_month), int(start_day), int(start_hour))
    end_dt = start_dt + timedelta(hours=Replace_Interval_Hours)
    end_year = str(end_dt.year)
    end_month = str(end_dt.month).zfill(2)   # 补零确保两位数
    end_day = str(end_dt.day).zfill(2)
    end_hour = str(end_dt.hour).zfill(2)

    # 计算运行时间
    run_days = str(Replace_Interval_Hours//24)
    run_hours = str(Replace_Interval_Hours%24)

    # 新的time_control配置内容
    new_time_control = f"""&time_control
    run_days                            = {run_days},
    run_hours                           = {run_hours},
    run_minutes                         = 0,
    run_seconds                         = 0,
    start_year                          = {start_year},
    start_month                         = {start_month},
    start_day                           = {start_day},
    start_hour                          = {start_hour}, 
    start_minute                        = 00, 
    start_second                        = 00, 
    end_year                            = {end_year},
    end_month                           = {end_month},
    end_day                             = {end_day},
    end_hour                            = {end_hour}, 
    end_minute                          = 00, 
    end_second                          = 00, 
    interval_seconds                    = 21600
    input_from_file                     = .true.,.true.,
    history_interval                    = 60, 
    frames_per_outfile                  = 1, 
    restart                             = .true.,
    restart_interval                    = {Replace_Interval_Hours * 60},  
    io_form_history                     = 2
    io_form_restart                     = 2
    io_form_input                       = 2
    io_form_boundary                    = 2
    """
    # 读取原始文件内容
    with open(namelist_input_file_path, 'r') as f:
        content = f.read()

    # 替换time_control部分（使用正则表达式精确匹配）
    pattern = r'(&time_control\s*[\s\S]*?)(?=/\s*\n|\Z)'
    updated_content = re.sub(
        pattern,
        new_time_control,
        content,
        count=1,
        flags=re.DOTALL
    )
    # 写入更新后的内容
    with open(namelist_input_file_path, 'w') as f:
        f.write(updated_content)

    print("namelist.input replace successful!")
    pass

def check_job_status(job_id):
    try:
        output = subprocess.check_output(["bjobs", job_id]).decode("utf-8")
        lines = output.splitlines()
        # 检查是否有至少两行
        if len(lines) >= 2:
            # **************************** 根据服务器命令进行修改 *******************************
            # 获取第二行的第三列（状态字段）
            status = lines[1].split()[2]
            # *********************************************************************************
            return status
        else:
            # 如果行数不足，返回 None
            return None
    except subprocess.CalledProcessError as e:
        # 如果命令执行失败，返回 None
        print(f"Error checking job status: {e}")
        return None
    
def wait_for_job_completion(job_id, wait_seconds, max_Run_Minutes):
    start_run = time.time()
    # **************************** 根据服务器命令进行修改 *******************************
    kill_command = f'scancel {job_id}'
    # *********************************************************************************
    while True:
        status = check_job_status(job_id) 
        run_Minutes = (time.time() - start_run)/60
        print(f"status: {status}, job_id: {job_id}")

        if run_Minutes > max_Run_Minutes: # 由于服务器故障，模式运行时间远远超出正常范围，kill掉这个作业
            print(f"Job {job_id} timed out after {max_Run_Minutes} minutes")
            subprocess.run(kill_command, shell=True, check=True)
            print(f"Job {job_id} has been killed due to timeout.")
            IF_job_complete = False
            break

        else:
            if status == "EXIT":
                print(f"JOB ID : {job_id} , Integrate Abnormally: {date_str}")
                IF_job_complete = False
                break

            elif status == "PEND":
                time.sleep(60)  # 秒，等待排队时间

            elif status == "RUN":
                time.sleep(wait_seconds)  # 等待运行

                try:
                    check_command = f"grep '{Success_Complate_Str}' {output_file}"
                    check_output = subprocess.check_output(check_command, shell=True).decode("utf-8")
                    if Success_Complate_Str in check_output:  # 检查是否已经成功运行完成
                        end_run = time.time()
                        print(f"JOB ID : {job_id} , Successfully integrate : {date_str} , Elapse Time: {end_run - start_run}Seconds")
                        IF_job_complete = True
                        break
                except subprocess.CalledProcessError:
                    pass                    

            else: # 状态信号为空，进一步检测：
                with open(f"slurm-{job_id}.out", 'r') as f:
                    content = f.read()
                    if Bad_Complate_Str in content:
                        print(f"JOB ID : {job_id} , Integrate Abnormally: {date_str}")
                        IF_job_complete = False
                        break
                    else:
                        try:
                            check_command = f"grep '{Success_Complate_Str}' {output_file}"
                            check_output = subprocess.check_output(check_command, shell=True).decode("utf-8")
                            if Success_Complate_Str in check_output:  # 检查是否输出文件是否包含成功运行完成的字符串
                                end_run = time.time()
                                print(f"JOB ID : {job_id} , Successfully integrate : {date_str} , Elapse Time: {end_run - start_run}Seconds")
                                IF_job_complete = True
                                break
                        except subprocess.CalledProcessError:
                            pass
    return IF_job_complete

# ---------------------------------- main --------------------------------------
start = time.time()

DOMAIN_PATH = '/public23/home/sc60844/dongyf/data/work/SMP/Domain'
RUN_PATH = '/public23/home/sc60844/dongyf/model/WRF/WRF-4.1.2/SMP/WRF-1/run'
SM_data = "CSMX"
DEPTH_ID =0
Replace_Interval_Hours = 48    #hour，土壤湿度循环替换周期
max_Run_Minutes = 120 # 单次模式运行最长时间（minutes）
wait_for_Run_seconds = 5 *60 # 刷新模式运行等待时间(seconds)
# wait_for_Run_seconds = 30  # 刷新模式运行等待时间(seconds): 模式运行时间的1/5或1/3。例如时间间隔2天，积分需要50min，可以设置10~20min的等待时间
Year = 2021
start_date_str = f"{Year}-082100"
end_date_str = f"{Year}-083118"
SM_PATH = f""
output_file = f"rsl.out.0000"
Success_Complate_Str = 'SUCCESS COMPLETE WRF'
Bad_Complate_Str = 'Exited with exit code'
# 生成时间序列
start_date = datetime.strptime(start_date_str, "%Y-%m%d%H")
end_date = datetime.strptime(end_date_str, "%Y-%m%d%H")
dates = pd.date_range(start=start_date, end=end_date, freq=f"{Replace_Interval_Hours}H")

# 读取静态数据
dst_lat, dst_lon, landmask, smcmax = read_geo_info(DOMAIN_PATH)

max_retries = 3  # 最大重试次数
abort_processing = False  # 全局终止标志

for date in dates:
    if abort_processing:
        print("ABORT: Processing terminated due to previous job failures")
        break  # 停止整个循环
    
    date_str = f'{date.strftime("%Y%m%d%H")}'
    print(f"Processing {date_str}")
    # --------------------- 数据准备步骤（只执行一次）---------------------
    # 将插值好的数据替换到wrfrst文件当中
    Replace_SM_data(RUN_PATH, SM_PATH, SM_data, date_str, dst_lon, dst_lat, landmask, smcmax) 
    # 修改namelist文件
    replace_namelist(RUN_PATH, date_str, Replace_Interval_Hours)
    # 切换到运行目录
    os.chdir(RUN_PATH)
    
    # 作业控制循环
    job_completed = False
    retry_count = 0
    SubmitCommand = "sbatch sub.sh"
    while not job_completed: # False
        # 提交作业（第一次或重试时执行）
        try:
            output = subprocess.check_output(SubmitCommand, shell=True).decode("utf-8")
            job_id = output.split()[-1]
            print(f"Submitted job {job_id} (attempt {retry_count+1}/{max_retries})")
        except Exception as e: # 提交失败，尝试重新提交
            print(f"Job submission failed: {str(e)}")
            retry_count += 1
            if retry_count >= max_retries:
                print(f"ERROR: Max submission retries reached for {date_str}")
                abort_processing = True  # 设置全局终止标志
                break
            sleep(60)
            continue
        
        # 等待作业完成
        IF_job_complete = wait_for_job_completion(job_id, wait_for_Run_seconds, max_Run_Minutes)
        
        if IF_job_complete:
            job_completed = True
            print(f"Job {job_id} completed successfully")
        else:
            retry_count += 1
            if retry_count >= max_retries:
                print(f"ERROR: Job {job_id} failed after {max_retries} attempts")
                abort_processing = True  # 设置全局终止标志
                break
            print(f"Job {job_id} failed. Retrying...")
            sleep(30)
    
    # 检查是否需要终止整个循环
    if abort_processing:
        print("CRITICAL: Terminating workflow due to job failure")
        break  # 终止整个循环
    
    # 继续执行下一个时次
    print(f"Completed {date_str}, proceeding to next time step")

# 最终状态报告
if abort_processing:
    print(f"Workflow ABORTED during {date_str}")
else:
    print(f"Successfully processed all {len(dates)} time steps")


print(f" End Date : {datetime.today()}")
end = time.time()
print(f"Elapse Time: {end - start}Seconds")    
    
    