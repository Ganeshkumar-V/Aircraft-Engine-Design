# Aircraft Engine Design Software
# (c) 2022 Ganeshkumar V

import tkinter as tk

import matplotlib
import numpy as np
from PIL import ImageTk, Image
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from dependencies.Turbofan import TurboFan_without_AB
from dependencies.Mixed_flow_Turbofan import mixed_flow_turbofan_without_afterburner as MFTF

matplotlib.use("TkAgg")


# Classes and Function
class Root(tk.Tk):
    def __init__(self, *args, **kwargs):
        super(Root, self).__init__()
        self.title("Aircraft Engine Designer")
        self.minsize(1300, 750)

class Graph:
    def __init__(self, frame, Title, Xlabel, Ylabel, Legend, x, *argv):
        # Plot Graph
        self.fig = Figure(figsize=(5, 5), dpi=100)
        a = self.fig.add_subplot(111)
        for arg in argv:
            a.plot(x, arg)
        a.set_xlabel(Xlabel)
        a.set_ylabel(Ylabel)
        a.set_title(Title)
        a.legend(Legend)
        self.canvas = FigureCanvasTkAgg(self.fig, frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        # Show Toolbar
        self.toolbar = NavigationToolbar2Tk(self.canvas, frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

class XGraph:
    def __init__(self, frame):
        self.fig = Figure(figsize=(5, 5), dpi=100)
        self.Frame = frame
        self.canvas = FigureCanvasTkAgg(self.fig, self.Frame)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.Frame)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def __call__(self, x, y, Label):
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(x, y, marker='x', label=Label)

    def wrap_it_up(self, Xlabel, Ylabel):
        self.ax.set_xlabel(Xlabel)
        self.ax.set_ylabel(Ylabel)
        self.ax.legend()
        self.canvas.draw()
        self.toolbar.update()

    def clear(self):
        self.fig.clf()

def validate_entry(inp):
    try:
        float(inp)
    except:
        return False
    return True

def draw_turbojet():
    Engine_Canvas.itemconfig(core_e, state='normal')
    Engine_Canvas.itemconfig(fan_e, state='hidden')
    Engine_Canvas.itemconfig(Ramjet_e, state='hidden')
    Engine_Canvas.itemconfig(Mixed_e, state='hidden')

def draw_turbofan():
    Engine_Canvas.itemconfig(fan_e, state='normal')
    Engine_Canvas.itemconfig(core_e, state='normal')
    Engine_Canvas.itemconfig(Ramjet_e, state='hidden')
    Engine_Canvas.itemconfig(Mixed_e, state='hidden')
    alpha_value_E.config(state='normal')
    alpha_range_min_E.config(state='normal')
    alpha_range_max_E.config(state='normal')
    alpha_range_Steps_E.config(state='normal')
    alpha_value.config(state='normal')
    alpha_range.config(state='normal')

def draw_ramjet():
    Engine_Canvas.itemconfig(fan_e, state='hidden')
    Engine_Canvas.itemconfig(core_e, state='hidden')
    Engine_Canvas.itemconfig(Ramjet_e, state='normal')
    Engine_Canvas.itemconfig(Mixed_e, state='hidden')

def draw_mixed_flow():
    Engine_Canvas.itemconfig(fan_e, state='hidden')
    Engine_Canvas.itemconfig(core_e, state='hidden')
    Engine_Canvas.itemconfig(Ramjet_e, state='hidden')
    Engine_Canvas.itemconfig(Mixed_e, state='normal')
    alpha_value_E.config(state='disabled')
    alpha_range_min_E.config(state='disabled')
    alpha_range_max_E.config(state='disabled')
    alpha_range_Steps_E.config(state='disabled')
    alpha_value.config(state='disabled')
    alpha_range.config(state='disabled')

def Enable_alpha_Range():
    if int(alpha_val.get()) == 1:  # Enable Alpha Range
        alpha_range_min_E.config(state='normal')  # disabled, normal or readonly
        alpha_range_max_E.config(state='normal')
        alpha_range_Steps_E.config(state='normal')
        alpha_value_E.config(state='disabled')
    else:  # Disable Alpha Range
        alpha_range_min_E.config(state='disabled')  # disabled, normal or readonly
        alpha_range_max_E.config(state='disabled')
        alpha_range_Steps_E.config(state='disabled')
        alpha_value_E.config(state='normal')

def Enable_Comp_R_Range():
    if int(Comp_R_val.get()) == 1:  # Enable Compression Ratio Range
        Comp_R_range_min_E.config(state='normal')  # disabled, normal or readonly
        Comp_R_range_max_E.config(state='normal')
        Comp_R_range_Steps_E.config(state='normal')
        Comp_R_value_E.config(state='disabled')
    else:  # Disable Compression Ratio Range
        Comp_R_range_min_E.config(state='disabled')  # disabled, normal or readonly
        Comp_R_range_max_E.config(state='disabled')
        Comp_R_range_Steps_E.config(state='disabled')
        Comp_R_value_E.config(state='normal')

def Enable_Fan_R_Range():
    if int(Fan_R_val.get()) == 1:  # Enable Fan Pressure Ratio Range
        Fan_R_range_min_E.config(state='normal')  # disabled, normal or readonly
        Fan_R_range_max_E.config(state='normal')
        Fan_R_range_Steps_E.config(state='normal')
        Fan_R_value_E.config(state='disabled')
    else:  # Disable Fan Pressure Ratio Range
        Fan_R_range_min_E.config(state='disabled')  # disabled, normal or readonly
        Fan_R_range_max_E.config(state='disabled')
        Fan_R_range_Steps_E.config(state='disabled')
        Fan_R_value_E.config(state='normal')

def ShowEngine(Engine):
    Engine_Canvas.itemconfig(Mixed_e, state='hidden')
    Engine_Canvas.itemconfig(Mixed_e_AB, state='hidden')
    Engine_Canvas.itemconfig(core_e, state='hidden')
    Engine_Canvas.itemconfig(core_e_AB, state='hidden')
    Engine_Canvas.itemconfig(fan_e, state='hidden')
    Engine_Canvas.itemconfig(fan_e_AB, state='hidden')
    Engine_Canvas.itemconfig(Ramjet_e, state='hidden')
    Engine_Canvas.itemconfig(Engine, state='normal')

def Enable_AB():
    if int(AB_val.get() == 1):
        # After burner
        AB_exit_E.config(state='normal')
        if int(Engine_val.get() == 1):          # Turbojet Engine
            ShowEngine(core_e_AB)
        elif int(Engine_val.get() == 2):        # Ramjet Engine
            ShowEngine(Ramjet_e)
        elif int(Engine_val.get() == 3):        # Mixed Flow Turbofan Engine
            ShowEngine(Mixed_e_AB)
            alpha_value_E.config(state='disabled')
            alpha_range_min_E.config(state='disabled')
            alpha_range_max_E.config(state='disabled')
            alpha_range_Steps_E.config(state='disabled')
            alpha_value.config(state='disabled')
            alpha_range.config(state='disabled')
        else:                                   # Turbofan Engine
            ShowEngine(fan_e_AB)
            alpha_value_E.config(state='normal')
            alpha_range_min_E.config(state='normal')
            alpha_range_max_E.config(state='normal')
            alpha_range_Steps_E.config(state='normal')
            alpha_value.config(state='normal')
            alpha_range.config(state='normal')
    else:
        # No After burner
        AB_exit_E.config(state='disabled')
        if int(Engine_val.get() == 1):  # Turbojet Engine
            ShowEngine(core_e)
        elif int(Engine_val.get() == 2):  # Ramjet Engine
            ShowEngine(Ramjet_e)
        elif int(Engine_val.get() == 3):  # Mixed Flow Turbofan Engine
            ShowEngine(Mixed_e)
            alpha_value_E.config(state='disabled')
            alpha_range_min_E.config(state='disabled')
            alpha_range_max_E.config(state='disabled')
            alpha_range_Steps_E.config(state='disabled')
            alpha_value.config(state='disabled')
            alpha_range.config(state='disabled')
        else:                                # Turbofan Engine
            ShowEngine(fan_e)
            Enable_alpha_Range()
            alpha_value.config(state='normal')
            alpha_range.config(state='normal')


def generate_array(Value, Max, Min, Steps, Val):
    Arr = []
    if Val.get() == 0:
        Arr.append(float(Value.get()))
    else:
        _Min_ = float(Min.get())
        _Max_ = float(Max.get())
        _Steps_ = int(Steps.get())
        Arr.append(_Min_)
        dx = (_Max_-_Min_)/_Steps_
        for i in range(1, _Steps_):
            Arr.append(_Min_ + dx * i)
        if Arr[len(Arr)-1] != _Max_:
            Arr.append(_Max_)
    nparr = np.array(Arr)
    return nparr

def solve():
    TL = int(Tlevel.get())
    M0 = float(Speed_E.get())
    h = float(Altitude_E.get())
    Fuel_name = Fuels.get()
    d = float(Diffuser_dia.get())
    Engine = int(Engine_val.get())
    AB = int(AB_val.get())
    if Engine == 0:
        bypass = generate_array(alpha_value_E, alpha_range_max_E, alpha_range_min_E, alpha_range_Steps_E, alpha_val)
    else:
        bypass = 0
    pi_c = generate_array(Comp_R_value_E, Comp_R_range_max_E, Comp_R_range_min_E, Comp_R_range_Steps_E, Comp_R_val)
    pi_f = generate_array(Fan_R_value_E, Fan_R_range_max_E, Fan_R_range_min_E, Fan_R_range_Steps_E, Fan_R_val)
    return pi_c, pi_f, TL, M0, h, Fuel_name, bypass, d

def enter():
    Welcome_Page.pack_forget()
    Input_Page.pack()
    # Analysis_Page.pack()
    # Input_Canvas.pack()
    # FL_Frame.pack()

def go_back():
    Analysis_Page.pack_forget()
    Input_Page.pack()

def generate_plots_comp(pi_c, pi_f, TL, M0, h, Fuel_name, bypass, d):
    length = len(pi_c)
    ST = np.zeros(length)
    TSFC = np.zeros(length)
    eta_p = np.zeros(length)
    eta_th = np.zeros(length)
    ST_G.clear()
    TSFC_G.clear()
    Eff_G.clear()
    for i in bypass:
        k = 0
        for j in pi_c:
            [Th, F, eta, Pt, Tt] = TurboFan_without_AB(j, pi_f, TL, M0, h, Fuel_name, i, d)
            ST[k] = Th[0]
            TSFC[k] = F[1]*3600
            eta_p[k] = eta[1]
            eta_th[k] = eta[0]
            k = k + 1

        # Plot Specific Thrust
        Label = "alpha: "+str(i)
        print(Label)
        ST_G(pi_c, ST, Label)
        TSFC_G(pi_c, TSFC, Label)
        Eff_G(pi_c, eta_p, Label)
    ST_G.wrap_it_up("Compressor Pressure Ratio", "Specific Thrust (N/kg/s)")
    TSFC_G.wrap_it_up("Compressor Pressure Ratio", "TSFC (kg/hr/N)")
    Eff_G.wrap_it_up("Compressor Pressure Ratio", "Propulsive Efficiency")

def generate_plots_MFTF(pi_c, pi_f, TL, M0, h, Fuel_name, d):
    length = len(pi_c)
    ST = np.zeros(length)
    TSFC = np.zeros(length)
    eta_p = np.zeros(length)
    eta_th = np.zeros(length)
    ST_G.clear()
    TSFC_G.clear()
    Eff_G.clear()
    for i in pi_f:
        k = 0
        for j in pi_c:
            [Th, F, eta, Pt, Tt] = MFTF(j, i, TL, M0, h, Fuel_name, d)
            ST[k] = Th[0]
            TSFC[k] = F[1]*3600
            eta_p[k] = eta[1]
            eta_th[k] = eta[0]
            k = k + 1

        # Plot Specific Thrust
        Label = "Fan Ratio: "+str(i)
        print(Label)
        ST_G(pi_c, ST, Label)
        TSFC_G(pi_c, TSFC, Label)
        Eff_G(pi_c, eta_p, Label)
    ST_G.wrap_it_up("Compressor Pressure Ratio", "Specific Thrust (N/kg/s)")
    TSFC_G.wrap_it_up("Compressor Pressure Ratio", "TSFC (kg/hr/N)")
    Eff_G.wrap_it_up("Compressor Pressure Ratio", "Propulsive Efficiency")

def generate_plots_fan(pi_c, pi_f, TL, M0, h, Fuel_name, bypass, d):
    length = len(pi_f)
    ST = np.zeros(length)
    TSFC = np.zeros(length)
    eta_p = np.zeros(length)
    eta_th = np.zeros(length)
    ST_G.clear()
    TSFC_G.clear()
    Eff_G.clear()
    for i in bypass:
        k = 0
        for j in pi_f:
            [Th, F, eta, Pt, Tt] = TurboFan_without_AB(pi_c, j, TL, M0, h, Fuel_name, i, d)
            ST[k] = Th[0]
            TSFC[k] = F[1]*3600
            eta_p[k] = eta[1]
            eta_th[k] = eta[0]
            k = k + 1

        # Plot Specific Thrust
        Label = "alpha: "+str(i)
        ST_G(pi_f, ST, Label)
        TSFC_G(pi_f, TSFC, Label)
        Eff_G(pi_f, eta_p, Label)
    ST_G.wrap_it_up("Fan Pressure Ratio", "Specific Thrust (N/kg/s)")
    TSFC_G.wrap_it_up("Fan Pressure Ratio", "TSFC (kg/hr/N)")
    Eff_G.wrap_it_up("Fan Pressure Ratio", "Propulsive Efficiency")

def ratio_selector(Pi_f_val, Pi_c_val, pi_c, pi_f, TL, M0, h, Fuel_name, bypass, d):
    if analysis_type.get() == 1:
        Ana_Comp_Frame.pack_forget()
        Ana_Fan_Frame.pack_forget()
        Graphs.pack_forget()
        Back.pack_forget()
        fan_ratio = float(Pi_f_val.get())
        print("Fan Pressure Ratio: ", fan_ratio)
        generate_plots_comp(pi_c, fan_ratio, TL, M0, h, Fuel_name, bypass, d)
        Ana_Comp_Frame.pack()
        Graphs.pack()
        Back.pack()
    elif analysis_type.get() == 2:
        Ana_Comp_Frame.pack_forget()
        Ana_Fan_Frame.pack_forget()
        Graphs.pack_forget()
        Back.pack_forget()
        comp_ratio = float(Pi_c_val.get())
        print("Compressor Pressure Ratio: ", comp_ratio)
        generate_plots_fan(comp_ratio, pi_f, TL, M0, h, Fuel_name, bypass, d)
        Ana_Fan_Frame.pack()
        Graphs.pack()
        Back.pack()
    return 0

def change_command():
    pi_c, pi_f, TL, M0, h, Fuel_name, alpha, d = solve()
    ratio_selector(Pi_f_val, Pi_c_val, pi_c, pi_f, TL, M0, h, Fuel_name, alpha, d)

def Effect_Compressor(pi_c, pi_f, TL, M0, h, Fuel_name, bypass, d):
    print("Inside Effect Compressor ", analysis_type.get())
    if int(analysis_type.get()) == 1:
        Ana_Comp_Frame.pack_forget()
        Ana_Fan_Frame.pack_forget()
        Back.pack_forget()
        Pi_f_options = list(map(str, pi_f))
        Pi_f_val.set(str(pi_f[0]))
        tk.Label(Ana_Comp_Frame, text="Fan Pressure Ratio : ").grid(row=0, column=0)
        tk.OptionMenu(Ana_Comp_Frame, Pi_f_val, *Pi_f_options,
                      command=ratio_selector(Pi_f_val, Pi_c_val, pi_c, pi_f, TL, M0, h, Fuel_name, bypass, d)).grid(row=0, column=1)
        tk.Button(Ana_Comp_Frame, text=" Change ", command=change_command).grid(row=0, column=2)
        Ana_Comp_Frame.pack()
    elif int(analysis_type.get()) == 2:
        Ana_Comp_Frame.pack_forget()
        Ana_Fan_Frame.pack_forget()
        Back.pack_forget()
        Pi_c_options = list(map(str, pi_c))
        Pi_c_val.set(Pi_c_options[0])
        print("Effect Compressor")
        tk.Label(Ana_Fan_Frame, text="Compressor Pressure Ratio : ").grid(row=0, column=0)
        tk.OptionMenu(Ana_Fan_Frame, Pi_c_val, *Pi_c_options,
                      command=ratio_selector(Pi_f_val, Pi_c_val, pi_c, pi_f, TL, M0, h, Fuel_name, bypass, d)).grid(row=0, column=1)
        tk.Button(Ana_Fan_Frame, text=" Change ", command=change_command).grid(row=0, column=2)
        Ana_Fan_Frame.pack()

def plot_command():
    pi_c, pi_f, TL, M0, h, Fuel_name, alpha, d = solve()
    Effect_Compressor(pi_c, pi_f, TL, M0, h, Fuel_name, alpha, d)

def analysis():
    if Engine_val.get() == 0:       # TurboFan Calculations
        Input_Page.pack_forget()
        Graphs.pack_forget()
        Back.pack_forget()
        Results.pack_forget()
        Analysis_Page.pack_forget()
        # Solvers
        pi_c, pi_f, TL, M0, h, Fuel_name, alpha, d = solve()
        print(" Pressure Ratio: ", pi_c)
        print(" Fan Ratio: ", pi_f)
        print(" Bypass Ratio: ", alpha)
        if (alpha_val.get() == 0) and (Comp_R_val.get() == 0) and (Fan_R_val.get() == 0):
            Th_para, fuel, Eff, Pt, Tt = TurboFan_without_AB(pi_c, pi_f, TL, M0, h, Fuel_name, alpha, d)
            ######################## Results Page #####################
            ####### Results Frame
            # Thrust frame
            Thrust_frame = tk.Frame(Results)
            tk.Label(Thrust_frame, text="Specific Thrust = ", font=("Times", 12)).grid(row=0, column=0, sticky=tk.E)
            tk.Label(Thrust_frame, text=str(round(Th_para[0], 4))+' N/kg/s', font=("Times", 12)).grid(row=0, column=1, sticky=tk.W)
            tk.Label(Thrust_frame, text="Mass flow rate of air = ", font=("Times", 12)).grid(row=1, column=0, sticky=tk.E)
            tk.Label(Thrust_frame, text=str(round(Th_para[1], 4))+' kg/s', font=("Times", 12)).grid(row=1, column=1, sticky=tk.W)
            tk.Label(Thrust_frame, text="Total Thrust = ", font=("Times", 12)).grid(row=2, column=0, sticky=tk.E)
            tk.Label(Thrust_frame, text=str(round(Th_para[2]/1000, 4))+' kN', font=("Times", 12)).grid(row=2, column=1, sticky=tk.W)
            Thrust_frame.grid(row=0, column=0, padx=25, pady=10)

            # Fuel Consumption Frame
            fuel_frame = tk.Frame(Results)
            tk.Label(fuel_frame, text="TSFC = ", font=("Times", 12)).grid(row=0, column=0, sticky=tk.E)
            tk.Label(fuel_frame, text=str(round(fuel[1]*3600, 4))+' kg/hr/N', font=("Times", 12)).grid(row=0, column=1, sticky=tk.W)
            tk.Label(fuel_frame, text="Fuel-Air Ratio = ", font=("Times", 12)).grid(row=1, column=0, sticky=tk.E)
            tk.Label(fuel_frame, text=str(round(fuel[0], 4)), font=("Times", 12)).grid(row=1, column=1, sticky=tk.W)
            tk.Label(fuel_frame, text="Bypass Ratio = ", font=("Times", 12)).grid(row=2, column=0, sticky=tk.E)
            tk.Label(fuel_frame, text=str(round(float(alpha_value_E.get()), 4)), font=("Times", 12)).grid(row=2, column=1, sticky=tk.W)
            fuel_frame.grid(row=0, column=1, padx=25, pady=10)

            # Efficiency Frame
            eff_frame = tk.Frame(Results)
            tk.Label(eff_frame, text="Propulsive Efficiency = ", font=("Times", 12)).grid(row=0, column=0, sticky=tk.E)
            tk.Label(eff_frame, text=str(round(Eff[1]*100, 4))+' %', font=("Times", 12)).grid(row=0, column=1, sticky=tk.W)
            tk.Label(eff_frame, text="Thermal Efficiency = ", font=("Times", 12)).grid(row=1, column=0, sticky=tk.E)
            tk.Label(eff_frame, text=str(round(Eff[0]*100, 4))+' %', font=("Times", 12)).grid(row=1, column=1, sticky=tk.W)
            tk.Label(eff_frame, text="Overall Efficiency = ", font=("Times", 12)).grid(row=2, column=0, sticky=tk.E)
            tk.Label(eff_frame, text=str(round(Eff[2]*100, 4))+' %', font=("Times", 12)).grid(row=2, column=1, sticky=tk.W)
            eff_frame.grid(row=0, column=2, padx=25, pady=10)
            Results.pack(pady=50)
            Back.pack(pady=20)
        else:
            ########################## Results Page ##########################
            ## Forgetting Previous Builds
            Results.pack_forget()
            Back.pack_forget()
            Graphs.pack_forget()

            ## Selection Analysis Frame
            Selection_Frame.pack()
            tk.Radiobutton(Selection_Frame, text="Effect of Compression Ratio",
                           variable=analysis_type, value=1).grid(row=0, column=0)
            tk.Radiobutton(Selection_Frame, text="Effect of Fan Pressure Ratio",
                           variable=analysis_type, value=2).grid(row=0, column=1)

        Analysis_Page.pack()
    elif Engine_val.get() == 3:    # Mixed Flow TurboFan
        # Forgetting Previous Frames
        Input_Page.pack_forget()
        Graphs.pack_forget()
        Back.pack_forget()
        Selection_Frame.pack_forget()
        Ana_Fan_Frame.pack_forget()
        Ana_Comp_Frame.pack_forget()
        Analysis_Page.pack_forget()
        # Get Inputs
        pi_c, pi_f, TL, M0, h, Fuel_name, unwanted, d = solve()
        if (Comp_R_val.get() == 0) and (Fan_R_val.get() == 0):
            Th_para, fuel, Eff, Pt, Tt = MFTF(pi_c, pi_f, TL, M0, h, Fuel_name, d)
            ######################## Results Page #####################
            ####### Results Frame
            # Thrust frame
            Thrust_frame = tk.Frame(Results)
            tk.Label(Thrust_frame, text="Specific Thrust = ", font=("Times", 12)).grid(row=0, column=0, sticky=tk.E)
            tk.Label(Thrust_frame, text=str(round(Th_para[0], 4)) + ' N/kg/s', font=("Times", 12)).grid(row=0, column=1,
                                                                                                      sticky=tk.W)
            tk.Label(Thrust_frame, text="Mass flow rate of air = ", font=("Times", 12)).grid(row=1, column=0,
                                                                                             sticky=tk.E)
            tk.Label(Thrust_frame, text=str(round(Th_para[1], 4)) + ' kg/s', font=("Times", 12)).grid(row=1, column=1,
                                                                                                      sticky=tk.W)
            tk.Label(Thrust_frame, text="Total Thrust = ", font=("Times", 12)).grid(row=2, column=0, sticky=tk.E)
            tk.Label(Thrust_frame, text=str(round(Th_para[2] / 1000, 4)) + ' kN', font=("Times", 12)).grid(row=2,
                                                                                                           column=1,
                                                                                                           sticky=tk.W)
            Thrust_frame.grid(row=0, column=0, padx=25, pady=10)

            # Fuel Consumption Frame
            fuel_frame = tk.Frame(Results)
            tk.Label(fuel_frame, text="TSFC = ", font=("Times", 12)).grid(row=0, column=0, sticky=tk.E)
            tk.Label(fuel_frame, text=str(round(fuel[1] * 3600, 4)) + ' kg/hr/N', font=("Times", 12)).grid(row=0,
                                                                                                           column=1,
                                                                                                           sticky=tk.W)
            tk.Label(fuel_frame, text="Fuel-Air Ratio = ", font=("Times", 12)).grid(row=1, column=0, sticky=tk.E)
            tk.Label(fuel_frame, text=str(round(fuel[0], 4)), font=("Times", 12)).grid(row=1, column=1, sticky=tk.W)
            tk.Label(fuel_frame, text="Bypass Ratio = ", font=("Times", 12)).grid(row=2, column=0, sticky=tk.E)
            tk.Label(fuel_frame, text=str(round(Th_para[3], 4)), font=("Times", 12)).grid(row=2, column=1, sticky=tk.W)
            fuel_frame.grid(row=0, column=1, padx=25, pady=10)

            # Efficiency Frame
            eff_frame = tk.Frame(Results)
            tk.Label(eff_frame, text="Propulsive Efficiency = ", font=("Times", 12)).grid(row=0, column=0, sticky=tk.E)
            tk.Label(eff_frame, text=str(round(Eff[1] * 100, 4)) + ' %', font=("Times", 12)).grid(row=0, column=1,
                                                                                                  sticky=tk.W)
            tk.Label(eff_frame, text="Thermal Efficiency = ", font=("Times", 12)).grid(row=1, column=0, sticky=tk.E)
            tk.Label(eff_frame, text=str(round(Eff[0] * 100, 4)) + ' %', font=("Times", 12)).grid(row=1, column=1,
                                                                                                  sticky=tk.W)
            tk.Label(eff_frame, text="Overall Efficiency = ", font=("Times", 12)).grid(row=2, column=0, sticky=tk.E)
            tk.Label(eff_frame, text=str(round(Eff[2] * 100, 4)) + ' %', font=("Times", 12)).grid(row=2, column=1,
                                                                                                  sticky=tk.W)
            eff_frame.grid(row=0, column=2, padx=25, pady=10)
            Results.pack(pady=50)
            Back.pack(pady=20)
        else:
            ############################## Results with Graph ###################
            ##### Forgetting Previous Builds
            Results.pack_forget()
            Back.pack_forget()
            Selection_Frame.pack_forget()
            Graphs.pack_forget()
            generate_plots_MFTF(pi_c, pi_f, TL, M0, h, Fuel_name, d)
            Graphs.pack()
            Back.pack(pady=20)
        Analysis_Page.pack()



# Main Program
root = Root()

# Variables
H = 12  # Altitude in km
# M0 = 2.5  # Flight Speed
# TL = 4  # Technology Level
Engine = 0  # 1 - Turbojet, 2 - Ramjet, 0 - TurboFan, 3 - Mixed flow Turbofan
AB = 0  # 0 - without After burner, 1 - with After Burner

################# Welcome Page #####################
# Title
Welcome_Page = tk.Frame(root, width=1300, height=700)
Title = tk.Label(
    Welcome_Page,
    text="Aircraft Engine Design",
    font=("Helvetica", 16, "bold italic"))  # font = ("Times", 16, "bold")
# Image
Welcome_Image = tk.Canvas(Welcome_Page, width=860, height=368)
img = ImageTk.PhotoImage(Image.open("images\Aircraft_image.png"))
Welcome_Image.create_image(0, 0, anchor=tk.NW, image=img)
# Enter button
Enter = tk.Button(Welcome_Page,
                  text="START ENGINE DESIGN",
                  command=enter,
                  bd=3,
                  width=20, height=1,
                  font=("Helvetica", 16, "bold"))

# Packing in the First Page
Title.grid(row=0, column=0, padx=0, pady=20)
Welcome_Image.grid(row=1, column=0, padx=0, pady=50)
Enter.grid(row=2, column=0, padx=0, pady=20)
Welcome_Page.pack()

###################### Input Page ######################
Input_Page = tk.Frame(root)
Fl_Tech_Frame = tk.Frame(Input_Page)
Engine_selection_Frame = tk.Frame(Input_Page)
Engine_Frame = tk.Frame(Input_Page)
Engine_Canvas = tk.Canvas(Input_Page)
Thermo_Frame = tk.Frame(Input_Page)
Perform_Analysis = tk.Frame(Input_Page)

# Title in Input Page
Title_in = tk.Label(
    Input_Page,
    text="Aircraft Engine Design",
    font=("Helvetica", 16, "bold italic"))  # font = ("Times", 16, "bold")
Title_in.pack(pady=10)

# Flight Conditions Frame
FL_Frame = tk.Frame(Fl_Tech_Frame, width=500)
Flight_conditions = tk.Label(FL_Frame, text="Flight Conditions ", font=("Helvetica", 11, "bold"))
Speed = tk.Label(FL_Frame, text="Flight Speed (M0) = ")
Speed_E = tk.Entry(FL_Frame, bd=2, validate="key", vcmd=(root.register(validate_entry), '%P'))
Altitude = tk.Label(FL_Frame, text="Altitude (km) = ")
Altitude_E = tk.Entry(FL_Frame, bd=2, validate="key", vcmd=(root.register(validate_entry), '%P'))
Flight_conditions.grid(row=0, column=0)
Speed.grid(row=1, column=0, sticky=tk.E)
Speed_E.grid(row=1, column=1, sticky=tk.W)
Altitude.grid(row=2, column=0, sticky=tk.E)
Altitude_E.grid(row=2, column=1, sticky=tk.W)
FL_Frame.grid(row=0, column=0, padx=50, pady=20)

# Technology Level Frame
Tech_Frame = tk.Frame(Fl_Tech_Frame, width=500)
Technology = tk.Label(Tech_Frame, text="Technology Level ", font=("Helvetica", 11, "bold"))
Technology.grid(row=0, column=0, sticky=tk.W)
Tlevel = tk.IntVar()
Tlevel.set(4)
T1 = tk.Radiobutton(Tech_Frame, text="Level 1", variable=Tlevel, value=1)
T2 = tk.Radiobutton(Tech_Frame, text="Level 2", variable=Tlevel, value=2)
T3 = tk.Radiobutton(Tech_Frame, text="Level 3", variable=Tlevel, value=3)
T4 = tk.Radiobutton(Tech_Frame, text="Level 4", variable=Tlevel, value=4)
T1.grid(row=1, column=0)
T2.grid(row=1, column=1)
T3.grid(row=2, column=0)
T4.grid(row=2, column=1)
Tech_Frame.grid(row=0, column=2, padx=50, pady=10)

# Engine Selection Frame
AB_F = tk.Frame(Engine_selection_Frame)
AB_label = tk.Label(AB_F, text="After burner ", font=("Helvetica", 11, "bold"))
AB_val = tk.IntVar()
AB1 = tk.Radiobutton(AB_F, text="without After Burner", variable=AB_val, value=0, command=Enable_AB)
AB2 = tk.Radiobutton(AB_F, text="with After Burner", variable=AB_val, value=1, command=Enable_AB)
Engine_F = tk.Frame(Fl_Tech_Frame)
Engine_label = tk.Label(Engine_F, text="Engine Configuration ", font=("Helvetica", 11, "bold"))
Engine_label.grid(row=0, column=0, sticky=tk.W)
Engine_val = tk.IntVar()
E1 = tk.Radiobutton(Engine_F, text="Turbojet", variable=Engine_val, value=1, command=Enable_AB)
E2 = tk.Radiobutton(Engine_F, text="Ramjet", variable=Engine_val, value=2, command=Enable_AB)
E3 = tk.Radiobutton(Engine_F, text="Turbofan", variable=Engine_val, value=0, command=Enable_AB)
E4 = tk.Radiobutton(Engine_F, text="Mixed flow Turbofan", variable=Engine_val, value=3, command=Enable_AB)
AB2.config(state='disabled')
E2.config(state='disabled')
E1.config(state='disabled')
E1.grid(row=1, column=0, sticky=tk.W)
E2.grid(row=1, column=1, sticky=tk.W)
E3.grid(row=2, column=0, sticky=tk.W)
E4.grid(row=2, column=1, sticky=tk.W)
Engine_F.grid(row=0, column=1, sticky=tk.W, padx=50)

Fl_Tech_Frame.pack(pady=5)

AB_label.grid(row=0, column=0, sticky=tk.W, padx=5)
AB1.grid(row=0, column=1, sticky=tk.W, padx=5)
AB2.grid(row=0, column=2, sticky=tk.W, padx=5)
AB_F.grid(row=0, column=0, padx=5, pady=5)
Engine_selection_Frame.pack()

# Input Canvas
# Turbojet Engine Schematic
core_engine = ImageTk.PhotoImage(Image.open("images\Turbojet_wo_AB.png"))
core_e = Engine_Canvas.create_image(0, 0, anchor=tk.NW, image=core_engine)
Engine_Canvas.itemconfig(core_e, state='hidden')
# Turbojet with AB Engine Schematic
Turbojet_engine_AB = ImageTk.PhotoImage(Image.open("images\Turbojet_w_AB.png"))
core_e_AB = Engine_Canvas.create_image(0, 0, anchor=tk.NW, image=Turbojet_engine_AB)
Engine_Canvas.itemconfig(core_e_AB, state='hidden')
# TurboFan Engine Schematic
fan = ImageTk.PhotoImage(Image.open("images\Turbofan_wo_AB.png"))
fan_e = Engine_Canvas.create_image(0, 0, anchor=tk.NW, image=fan)
Engine_Canvas.itemconfig(fan_e, state='normal')
# Turbofan with AB Schematic
fan_AB = ImageTk.PhotoImage(Image.open("images\Turbofan_w_AB.png"))
fan_e_AB = Engine_Canvas.create_image(0, 0, anchor=tk.NW, image=fan_AB)
Engine_Canvas.itemconfig(fan_e_AB, state='hidden')
# Mixed TurboFan Engine Schematic
Mixed_engine = ImageTk.PhotoImage(Image.open("images\Mixed_flow_wo_AB.png"))
Mixed_e = Engine_Canvas.create_image(0, 0, anchor=tk.NW, image=Mixed_engine)
Engine_Canvas.itemconfig(Mixed_e, state='hidden')
# Mixed Flow TurboFan with AfterBurner Engine Schematic
Mixed_engine_AB = ImageTk.PhotoImage(Image.open("images\Mixed_flow_w_AB.png"))
Mixed_e_AB = Engine_Canvas.create_image(0, 0, anchor=tk.NW, image=Mixed_engine_AB)
Engine_Canvas.itemconfig(Mixed_e_AB, state='hidden')
# Ramjet Engine Schematic
Ramjet_engine = ImageTk.PhotoImage(Image.open("images\Ramjet.png"))
Ramjet_e = Engine_Canvas.create_image(0, 0, anchor=tk.NW, image=Ramjet_engine)
Engine_Canvas.itemconfig(Ramjet_e, state='hidden')
Engine_Canvas.pack(padx=0, pady=20)

# Engine Parameters
# Bypass Ratio Inputs
TurboFan = tk.Frame(Engine_Frame)
alpha_F = tk.Frame(TurboFan)
alpha_val = tk.IntVar()
alpha_label = tk.Label(alpha_F, text="Bypass Ratio: ", font=("Helvetica", 10, "bold"))
alpha_label.grid(row=0, column=0)
alpha_range_F = tk.Frame(alpha_F)
alpha_range_min = tk.Label(alpha_range_F, text="Min ")
alpha_range_min.grid(row=0, column=0)
alpha_range_min_E = tk.Entry(alpha_range_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
alpha_range_min_E.grid(row=0, column=1)
alpha_range_max = tk.Label(alpha_range_F, text="Max ")
alpha_range_max.grid(row=1, column=0)
alpha_range_max_E = tk.Entry(alpha_range_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
alpha_range_max_E.grid(row=1, column=1)
alpha_range_Steps = tk.Label(alpha_range_F, text="Steps ")
alpha_range_Steps.grid(row=2, column=0)
alpha_range_Steps_E = tk.Entry(alpha_range_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
alpha_range_Steps_E.grid(row=2, column=1)
alpha_range_F.grid(row=1, column=2)
alpha_value_E = tk.Entry(alpha_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
alpha_value_E.grid(row=0, column=2)
alpha_value = tk.Radiobutton(alpha_F, text="Value", variable=alpha_val, value=0, command=Enable_alpha_Range)
alpha_value.grid(row=0, column=1)
alpha_range = tk.Radiobutton(alpha_F, text="Range", variable=alpha_val, value=1, command=Enable_alpha_Range)
alpha_range_min_E.config(state='disabled')
alpha_range_max_E.config(state='disabled')
alpha_range_Steps_E.config(state='disabled')
alpha_value_E.config(state='normal')
alpha_range.grid(row=1, column=1)
alpha_F.grid(row=0, column=0, padx=10, pady=5)

# Compressor Pressure Ratio Inputs
Comp_RF = tk.Frame(TurboFan)
Comp_R_val = tk.IntVar()
Comp_R_label = tk.Label(Comp_RF, text="Compressor Ratio: ", font=("Helvetica", 10, "bold"))
Comp_R_label.grid(row=0, column=0)
Comp_R_range_F = tk.Frame(Comp_RF)
Comp_R_range_min = tk.Label(Comp_R_range_F, text="Min ")
Comp_R_range_min.grid(row=0, column=0)
Comp_R_range_min_E = tk.Entry(Comp_R_range_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
Comp_R_range_min_E.grid(row=0, column=1)
Comp_R_range_max = tk.Label(Comp_R_range_F, text="Max ")
Comp_R_range_max.grid(row=1, column=0)
Comp_R_range_max_E = tk.Entry(Comp_R_range_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
Comp_R_range_max_E.grid(row=1, column=1)
Comp_R_range_Steps = tk.Label(Comp_R_range_F, text="Steps ")
Comp_R_range_Steps.grid(row=2, column=0)
Comp_R_range_Steps_E = tk.Entry(Comp_R_range_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
Comp_R_range_Steps_E.grid(row=2, column=1)
Comp_R_range_F.grid(row=1, column=2)
Comp_R_value = tk.Radiobutton(Comp_RF, text="Value", variable=Comp_R_val, value=0, command=Enable_Comp_R_Range)
Comp_R_value_E = tk.Entry(Comp_RF, validate="key", vcmd=(root.register(validate_entry), '%P'))
Comp_R_value_E.grid(row=0, column=2)
Comp_R_range = tk.Radiobutton(Comp_RF, text="Range", variable=Comp_R_val, value=1, command=Enable_Comp_R_Range)
Comp_R_range_min_E.config(state='disabled')  # disabled, normal or readonly
Comp_R_range_max_E.config(state='disabled')
Comp_R_range_Steps_E.config(state='disabled')
Comp_R_value.grid(row=0, column=1)
Comp_R_range.grid(row=1, column=1)
Comp_RF.grid(row=0, column=1, padx=10, pady=5)

# Fan Pressure Ratio Inputs
Fan_RF = tk.Frame(TurboFan)
Fan_R_val = tk.IntVar()
Fan_R_label = tk.Label(Fan_RF, text="Fan Pressure Ratio: ", font=("Helvetica", 10, "bold"))
Fan_R_label.grid(row=0, column=0)
Fan_R_range_F = tk.Frame(Fan_RF)
Fan_R_range_min = tk.Label(Fan_R_range_F, text="Min ")
Fan_R_range_min.grid(row=0, column=0)
Fan_R_range_min_E = tk.Entry(Fan_R_range_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
Fan_R_range_min_E.grid(row=0, column=1)
Fan_R_range_max = tk.Label(Fan_R_range_F, text="Max ")
Fan_R_range_max.grid(row=1, column=0)
Fan_R_range_max_E = tk.Entry(Fan_R_range_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
Fan_R_range_max_E.grid(row=1, column=1)
Fan_R_range_Steps = tk.Label(Fan_R_range_F, text="Steps ")
Fan_R_range_Steps.grid(row=2, column=0)
Fan_R_range_Steps_E = tk.Entry(Fan_R_range_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
Fan_R_range_Steps_E.grid(row=2, column=1)
Fan_R_range_F.grid(row=1, column=2)
Fan_R_value = tk.Radiobutton(Fan_RF, text="Value", variable=Fan_R_val, value=0, command=Enable_Fan_R_Range)
Fan_R_value_E = tk.Entry(Fan_RF, validate="key", vcmd=(root.register(validate_entry), '%P'))

Fan_R_value_E.grid(row=0, column=2)
Fan_R_range = tk.Radiobutton(Fan_RF, text="Range", variable=Fan_R_val, value=1, command=Enable_Fan_R_Range)
Fan_R_range_min_E.config(state='disabled')
Fan_R_range_max_E.config(state='disabled')
Fan_R_range_Steps_E.config(state='disabled')
Fan_R_value.grid(row=0, column=1, sticky=tk.W)
Fan_R_range.grid(row=1, column=1, sticky=tk.W)
Fan_RF.grid(row=0, column=2, padx=10, pady=5)
TurboFan.pack()

# Jet Fuel Frame
Jet_Fuel_F = tk.Frame(Thermo_Frame)
Fuel_Label = tk.Label(Jet_Fuel_F, text="Jet Fuel", font=("Helvetica", 10, "bold"))
Fuel_Options = ["JET A", "JET A1", "AVGAS 100L"]
Fuels = tk.StringVar()
Fuels.set(Fuel_Options[0])
Fuel_choose = tk.OptionMenu(Jet_Fuel_F, Fuels, *Fuel_Options)
Fuel_choose.config(width=15)
Fuel_Label.grid(row=0, column=0)
Fuel_choose.grid(row=0, column=1)
Jet_Fuel_F.grid(row=0, column=0, padx=30, pady=0)
# Turbine Inlet Temperature Frame
Temperature_F = tk.Frame(Thermo_Frame)
tk.Label(Temperature_F, text="Temperature (K)", font=("Helvetica", 10, "bold")).grid(row=0, column=0, sticky=tk.W)
tk.Label(Temperature_F, text="Turbine Inlet").grid(row=1, column=0, sticky=tk.E)
Turbine_E = tk.Entry(Temperature_F, validate="key", vcmd=(root.register(validate_entry), '%P')).grid(row=1, column=1,
                                                                                                     sticky=tk.W)
tk.Label(Temperature_F, text="AB Exit").grid(row=2, column=0, sticky=tk.E)
AB_exit_E = tk.Entry(Temperature_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
AB_exit_E.grid(row=2, column=1, sticky=tk.W)
AB_exit_E.config(state='disabled')
Temperature_F.grid(row=0, column=1, padx=30, pady=0)
# Diffuser Area
Diffuser_F = tk.Frame(Thermo_Frame)
tk.Label(Diffuser_F, text="Engine Size", font=("Helvetica", 10, "bold")).grid(row=0, column=0, sticky=tk.W)
tk.Label(Diffuser_F, text="Diffuser Inlet Diameter (m) ").grid(row=1, column=0, sticky=tk.W)
Diffuser_dia = tk.Entry(Diffuser_F, validate="key", vcmd=(root.register(validate_entry), '%P'))
Diffuser_dia.grid(row=1, column=1, sticky=tk.W)
Diffuser_F.grid(row=0, column=2, padx=30, pady=10)

Thermo_Frame.pack()
Engine_Frame.pack(pady=20)

# Testing Purposes
# Diffuser_dia.insert(0, "1.95")
# Fan_R_value_E.insert(0, "2")
# Speed_E.insert(0, "0.89")
# Altitude_E.insert(0, "12")
# alpha_value_E.insert(0, "2")
# Comp_R_value_E.insert(0, "20")

# Perform Analysis Frame
Analysis = tk.Button(Perform_Analysis,
                     text="PERFORM ANALYSIS",
                     command=analysis,
                     bd=3,
                     font=("Helvetica", 12, "bold"))
Analysis.grid(row=0, column=0, padx=0, pady=10)
Perform_Analysis.pack()
# print("Flight Speed is ", M0, " and Altitude is ", H, " Technology Level ", TL, " is selected!")

############################## Analysis page ###############################
Analysis_Page = tk.Frame(root)

####### Title in Analysis Page
Title_an = tk.Label(
    Analysis_Page,
    text="Aircraft Engine Analysis Results",
    font=("Helvetica", 16, "bold italic"))  # font = ("Times", 16, "bold")
Title_an.pack(pady=20)
## Effect of Compressor Pressure Ratio
analysis_type = tk.IntVar()
Pi_f_val = tk.StringVar()
Pi_c_val = tk.StringVar()
Selection_Frame = tk.Frame(Analysis_Page)
plot = tk.Button(Selection_Frame, text=" PLOT ", command=plot_command)
plot.grid(row=0, column=3)
Results = tk.Frame(Analysis_Page)
# Graph Frames
Graphs = tk.Frame(Analysis_Page)
ST_F = tk.Frame(Graphs)
TSFC_F = tk.Frame(Graphs)
Eff_F = tk.Frame(Graphs)
ST_G = XGraph(ST_F)
TSFC_G = XGraph(TSFC_F)
Eff_G = XGraph(Eff_F)
ST_F.grid(row=0, column=0, padx=10, pady=10)
TSFC_F.grid(row=0, column=1, padx=10, pady=10)
Eff_F.grid(row=0, column=2, padx=10, pady=10)
Ana_Comp_Frame = tk.Frame(Analysis_Page)
Ana_Fan_Frame = tk.Frame(Analysis_Page)

Back = tk.Button(Analysis_Page, text="GO BACK", command=go_back, bd=3, font=("Helvetica", 12, "bold"))
# Start of Main Loop
root.mainloop()

# Model Syntax

# Checking PNG
# canvas2 = Canvas(root, width = 600, height = 400)
# canvas2.pack()
# img2 = ImageTk.PhotoImage(Image.open("graphics.png"))
# canvas2.create_image(0, 0, anchor=NW, image=img2)

# Checking Arrow mark
# canvas2 = Canvas(root)
# canvas2.pack()
# canvas2.create_line(0, 0, 200, 100, arrow=LAST)
# canvas2.create_line(200, 100, 100, 0, arrow=LAST)
