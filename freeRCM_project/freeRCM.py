#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 12:59:02 2024

@author: Piero Wemyss

TO DO:

    -- Add help menus
    -- --DONE-ish--     Add click to gen line (only supports single lines at a time at the moment)
    -- Add Pitzer correlation acentricity generation
    -- get nuitka startup to go faster for standalone exec/.exe/app bundle
    -- Get rid of bloat/scratch code

"""

import sys
import os
import numpy as np
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QPushButton, QVBoxLayout, QHBoxLayout,
                               QLabel, QListWidget, QLineEdit, QComboBox, QFrame, QTableWidget,
                               QTableWidgetItem, QFileDialog, QMessageBox, QInputDialog, QGridLayout,
                               QSizePolicy, QSpacerItem)
from PySide6.QtCore import Qt, QSize
from PySide6.QtGui import QPixmap
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
import pickle
from dict2struct import dict2struct
from RCM import RCM
from RCMplot import RCMplot


global P, comps, selected_comps, allProps, lmopts, opts, NRTL_aij, NRTL_bij, NRTL_cij, TcCel, Pc, omega, antoine_params, PLXANT_params

P = 1
comps = np.array([])
selected_comps = np.array([])
antoine_params = np.array([])
PLXANT_params = np.array([])
NRTL_aij = np.array([])
NRTL_bij = np.array([])
NRTL_cij = np.array([])
TcCel = np.array([])
Pc = np.array([])
omega = np.array([])

lmopts = {
    'maxiter': 500,
    'ftol': 1e-12,
    'xtol': 1e-12
}
options = {
    "antMethod": 2,
    "activity": 3,
    "lines": 15,
    "linewidth": 1.2,
    "n_it": 80,
    "lmopts": lmopts
}
opts = dict2struct(options)

class GetStartedWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Get Started")
        self.setGeometry(100, 100, 600, 400)
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.create_widgets()

    def create_widgets(self):
        layout = QGridLayout(self.central_widget)

        welcome_label = QLabel("Welcome to FreeRCM", self)
        welcome_label.setStyleSheet("font-size: 30px;")
        layout.addWidget(welcome_label, 0, 0, 1, 2)
        
        def resource_path(relative_path):
            try:
                base_path = sys._MEIPASS
            except Exception:
                base_path = os.path.abspath(".")

            return os.path.join(base_path, relative_path)

        logo_path = resource_path('logo.png')

        logo_label = QLabel(self)
        logo_label.setPixmap(QPixmap(logo_path).scaled(500, 500, Qt.AspectRatioMode.KeepAspectRatio))
        layout.addWidget(logo_label, 0, 2, 3, 1)

        new_simulation_button = QPushButton("New Simulation", self)
        new_simulation_button.clicked.connect(self.open_new_simulation)
        layout.addWidget(new_simulation_button, 1, 0)

        open_simulation_button = QPushButton("Open Simulation", self)
        open_simulation_button.clicked.connect(self.open_simulation)
        layout.addWidget(open_simulation_button, 1, 1)

    def open_new_simulation(self):
        global P, comps, selected_comps, allProps, lmopts, opts, NRTL_aij, NRTL_bij, NRTL_cij, TcCel, Pc, omega, antoine_params, PLXANT_params
        P = 1
        comps = np.array([])
        selected_comps = np.array([])
        antoine_params = np.array([])
        PLXANT_params = np.array([])
        NRTL_aij = np.array([])
        NRTL_bij = np.array([])
        NRTL_cij = np.array([])
        TcCel = np.array([])
        Pc = np.array([])
        omega = np.array([])

        self.hide()
        self.new_simulation_window = NewSimulationWindow()
        self.new_simulation_window.show()

    def open_simulation(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select file", "", "RCM Files (*.rcm)")
        if file_path:
            global P, comps, selected_comps, allProps, lmopts, opts, NRTL_aij, NRTL_bij, NRTL_cij, TcCel, Pc, omega, antoine_params, PLXANT_params
            with open(file_path, 'rb') as file:
                data = pickle.load(file)
            P = data['P']
            comps = data['comps']
            selected_comps = data['selected_comps']
            allProps = data['allProps']
            lmopts = data['lmopts']
            opts = data['opts']
            NRTL_aij = data['NRTL_aij']
            NRTL_bij = data['NRTL_bij']
            NRTL_cij = data['NRTL_cij']
            TcCel = data['TcCel']
            Pc = data['Pc']
            omega = data['omega']
            antoine_params = data['antoine_params']
            PLXANT_params = data['PLXANT_params']
#             print(f"Loaded variables: P = {P}, comps = {comps}, ..., PLXANT_params = {PLXANT_params}")
            self.hide()
            self.new_simulation_window = NewSimulationWindow()
            self.new_simulation_window.show()


class NewSimulationWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        global comps, selected_comps
        if comps.size:
            self.setWindowTitle("Existing Simulation")
        else:
            self.setWindowTitle("New Simulation")
        self.setGeometry(100, 100, 800, 600)
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.create_widgets()

    def create_widgets(self):
        global opts
        layout = QGridLayout(self.central_widget)

        add_components_button = QPushButton("Add Components", self)
        add_components_button.clicked.connect(self.add_components)
        layout.addWidget(add_components_button, 0, 0)

        delete_components_button = QPushButton("Delete Components", self)
        delete_components_button.clicked.connect(self.delete_components)
        layout.addWidget(delete_components_button, 0, 1)

        delete_components_button.setMaximumWidth(add_components_button.sizeHint().width())

        midsection = QFrame()
        midsection_cont = QHBoxLayout(midsection)
        midsection_cont.setSpacing(10)
        
        left_column = QVBoxLayout()
        self.components_list_label = QLabel("All Components", self)
        self.components_list_label.setFixedHeight(20)
        left_column.addWidget(self.components_list_label)
        
        self.components_list = QListWidget(self)
        left_column.addWidget(self.components_list)
        
        midsection_cont.addLayout(left_column, 2)
        
        middle_column = QVBoxLayout()
        middle_column.addItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))
        move_to_selected_button = QPushButton(">>", self)
        move_to_selected_button.clicked.connect(self.move_to_selected)
        middle_column.addWidget(move_to_selected_button)
        
        move_to_master_button = QPushButton("<<", self)
        move_to_master_button.clicked.connect(self.move_to_master)
        middle_column.addWidget(move_to_master_button)
        middle_column.addItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))
        
        midsection_cont.addLayout(middle_column)
        
        right_column = QVBoxLayout()
        self.selected_comps_label = QLabel("Selected Components", self)
        self.selected_comps_label.setFixedHeight(20)
        right_column.addWidget(self.selected_comps_label)
        
        self.selected_comps_list = QListWidget(self)
        self.update_selected_components_list()
        self.update_components_list()
        right_column.addWidget(self.selected_comps_list)
        
        midsection_cont.addLayout(right_column, 2)
        
        layout.addWidget(midsection, 1, 0, 2, 5)

        help_button = QPushButton("Help", self)
        help_button.clicked.connect(self.show_help)
        layout.addWidget(help_button, 0, 4)

        vapor_pressure_label = QLabel("Vapor Pressure:", self)
        layout.addWidget(vapor_pressure_label, 3, 0)

        self.vapor_pressure_var = QComboBox(self)
        self.vapor_pressure_var.addItems(["Antoine", "Extended Antoine"])
        self.vapor_pressure_var.setCurrentIndex(opts.antMethod-1)
        self.vapor_pressure_var.currentIndexChanged.connect(self.update_vapor_pressure_method)
        layout.addWidget(self.vapor_pressure_var, 3, 1)

        dropdown_label = QLabel("Select Model:", self)
        layout.addWidget(dropdown_label, 4, 0)

        self.model_var = QComboBox(self)
        self.model_var.addItems(["Ideal", "NRTL", "NRTL-SRK"])
        self.model_var.setCurrentIndex(opts.activity-1)
        self.model_var.currentIndexChanged.connect(self.update_eos_method)
        layout.addWidget(self.model_var, 4, 1)

        self.input_parameters_button = QPushButton("Input Parameters", self)
        self.input_parameters_button.clicked.connect(self.input_parameters)
        layout.addWidget(self.input_parameters_button, 4, 2)

        pressure_label = QLabel("Pressure:", self)
        layout.addWidget(pressure_label, 5, 0)

        self.pressure_entry = QLineEdit(str(P), self)
        self.pressure_entry.returnPressed.connect(self.save_pressure)
        layout.addWidget(self.pressure_entry, 5, 1)

        pressure_unit_label = QLabel(" bar", self)
        layout.addWidget(pressure_unit_label, 5, 2)

        go_back_button = QPushButton("Go Back", self)
        go_back_button.clicked.connect(self.go_back)
        layout.addWidget(go_back_button, 6, 0)

        next_button = QPushButton("Next", self)
        next_button.clicked.connect(self.next)
        layout.addWidget(next_button, 6, 4)

        self.central_widget.setFocusPolicy(Qt.FocusPolicy.StrongFocus)


    def add_components(self):
        new_comp, ok = QInputDialog.getText(self, "Add Component", "Enter component name:")
        if ok and new_comp:
            global comps
            comps = np.append(comps, new_comp)
            self.update_components_list()

    def delete_components(self):
        selected_items = self.components_list.selectedItems()
        if selected_items:
            selected_item = selected_items[0]
            selected_index = self.components_list.row(selected_item)
            global comps
            comps = np.delete(comps, selected_index)
            self.update_components_list()

    def update_components_list(self):
        self.components_list.clear()
        self.components_list.addItems(comps.tolist())
        self.update_selected_components_list()

    def update_selected_components_list(self):
        global selected_comps
        self.selected_comps_list.clear()
        self.selected_comps_list.addItems(selected_comps.tolist())

    def move_to_selected(self):
        selected_items = self.components_list.selectedItems()
        if selected_items:
            global selected_comps, comps
            selected_item = selected_items[0]
            selected_comps = np.append(selected_comps, selected_item.text())
            self.update_components_list()

    def move_to_master(self):
        selected_items = self.selected_comps_list.selectedItems()
        if selected_items:
            selected_item = selected_items[0]
            global selected_comps, comps
            selected_comps = np.delete(selected_comps, self.selected_comps_list.row(selected_item))
            self.update_components_list()

    def save_pressure(self):
        global P
        try:
            P = float(self.pressure_entry.text())
        except ValueError:
            QMessageBox.critical(self, "Invalid Input", "Please enter a valid number for pressure.")

    def input_parameters(self):
        self.input_params_window = InputParamsWindow()
        self.input_params_window.show()

    def update_vapor_pressure_method(self):
        global vapor_pressure_method, opts
        vapor_pressure_method = self.vapor_pressure_var.currentText()
        if vapor_pressure_method == "Antoine":
            opts.antMethod = 1
        elif vapor_pressure_method == "Extended Antoine":
            opts.antMethod = 2

    def update_eos_method(self):
        global eos_method, opts
        eos_method = self.model_var.currentText()
        if eos_method == "Ideal":
            opts.activity = 1
        elif eos_method == "NRTL":
            opts.activity = 2
        elif eos_method == "NRTL-SRK":
            opts.activity = 3

    def show_help(self):
        QMessageBox.information(self, "Help", "This is a simulation tool.")

    def go_back(self):
        self.hide()
        self.get_started_window = GetStartedWindow()
        self.get_started_window.show()

    def next(self):
#         self.hide()
#         self.make_sim_window = MakeSimWindow()
#         self.make_sim_window.show()

        global selected_comps, opts, NRTL_aij, NRTL_bij, NRTL_cij, crit_params, antoine_params, PLXANT_params
        if len(selected_comps) == 3:
            if opts.antMethod == 1 and antoine_params.size:
                if opts.activity == 1: 
                    self.hide()
                    self.make_sim_window = MakeSimWindow()
                    self.make_sim_window.show()
                elif opts.activity == 2 and NRTL_aij.size and NRTL_bij.size and NRTL_cij.size:
                    self.hide()
                    self.make_sim_window = MakeSimWindow()
                    self.make_sim_window.show()
                elif opts.activity == 3 and NRTL_aij.size and NRTL_bij.size and NRTL_cij.size and TcCel.size and Pc.size and omega.size:
                    self.hide()
                    self.make_sim_window = MakeSimWindow()
                    self.make_sim_window.show()
                else:
                    QMessageBox.information(self, "Help", "Please ensure all parameters are inserted correctly for your selected fluid package.")
            elif opts.antMethod == 2 and PLXANT_params.size:
                if opts.activity == 1: 
                    self.hide()
                    self.make_sim_window = MakeSimWindow()
                    self.make_sim_window.show()
                elif opts.activity == 2 and NRTL_aij.size and NRTL_bij.size and NRTL_cij.size:
                    self.hide()
                    self.make_sim_window = MakeSimWindow()
                    self.make_sim_window.show()
                elif opts.activity == 3 and NRTL_aij.size and NRTL_bij.size and NRTL_cij.size and TcCel.size and Pc.size and omega.size:
                    self.hide()
                    self.make_sim_window = MakeSimWindow()
                    self.make_sim_window.show()
                else:
                    QMessageBox.information(self, "Help", "Please ensure all parameters are inserted correctly for your selected fluid package.")
            else:
                QMessageBox.information(self, "Help", "Please ensure all parameters are inserted correctly for your selected vapor pressure model")
        else:
            QMessageBox.information(self, "Help", "Please select only 3 components.")


class MakeSimWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        global P, comps, selected_comps, opts
        self.comps = comps
        self.selected_comps = selected_comps
        self.opts = opts
        self.setWindowTitle("Simulation - Auto-generate or click on plot!")
        self.create_widgets()
        
    def create_widgets(self):
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        
        self.main_layout = QHBoxLayout(self.central_widget)
        
        self.button_frame = QFrame()
        self.button_layout = QVBoxLayout(self.button_frame)
        
        self.auto_gen_button = QPushButton("Auto-Generate Curves")
        self.auto_gen_button.clicked.connect(self.auto_generate)
        self.button_layout.addWidget(self.auto_gen_button)
        
        self.clear_button = QPushButton("Clear Plot")
        self.clear_button.clicked.connect(self.clear_plot)
        self.button_layout.addWidget(self.clear_button)

        self.adj_params_button = QPushButton("Adjust Parameters")
        self.adj_params_button.clicked.connect(self.adjust_parameters)
        self.button_layout.addWidget(self.adj_params_button)
        
        self.plot_options_button = QPushButton("Plot Options")
        self.plot_options_button.clicked.connect(self.plot_options)
        self.button_layout.addWidget(self.plot_options_button)
        
        self.solver_options_button = QPushButton("Solver Options")
        self.solver_options_button.clicked.connect(self.solver_options)
        self.button_layout.addWidget(self.solver_options_button)
        
        self.save_button = QPushButton("Save Simulation")
        self.save_button.clicked.connect(self.save_variables)
        self.button_layout.addWidget(self.save_button)
        
        self.go_back_button = QPushButton("Go Back")
        self.go_back_button.clicked.connect(self.go_back)
        self.button_layout.addWidget(self.go_back_button)
        
        self.main_layout.addWidget(self.button_frame)
        
        self.plot_frame = QFrame()
        self.main_layout.addWidget(self.plot_frame)
        
        self.canvas = None
        self.clear_figure()
        # self.plot_figure()

    def plot_figure(self):
        global comps, selected_comps, P, allProps, opts
        x0n = np.array([])
        x = RCM(comps,selected_comps,P,allProps,opts,x0n,1).x
        fig, ax = RCMplot(x, self.selected_comps, self.opts)
        
        if self.canvas:
            self.canvas.deleteLater()
            self.toolbar.deleteLater()
        
        self.canvas = FigureCanvas(fig)
        self.toolbar = CompactNavigationToolbar(self.canvas, self)
        self.coord_label = QLabel("", self.toolbar)
        self.toolbar.addWidget(self.coord_label)
        
        self.plot_layout = QVBoxLayout(self.plot_frame)
        self.plot_layout.addWidget(self.toolbar)
        self.plot_layout.addWidget(self.canvas)
        
        self.canvas.draw()

        self.canvas.mpl_connect('button_press_event', self.click_plot)
 
    def clear_figure(self):
        x = np.zeros([1,3,1])
        fig, ax = RCMplot(x, self.selected_comps, self.opts)
        
        if self.canvas:
            self.canvas.deleteLater()
            self.toolbar.deleteLater()
        
        self.canvas = FigureCanvas(fig)
        self.toolbar = CompactNavigationToolbar(self.canvas, self)
        self.coord_label = QLabel("", self.toolbar)
        self.toolbar.addWidget(self.coord_label)
        
        self.plot_layout = QVBoxLayout(self.plot_frame)
        self.plot_layout.addWidget(self.toolbar)
        self.plot_layout.addWidget(self.canvas)
        
        self.canvas.draw()

        self.canvas.mpl_connect('button_press_event', self.click_plot)

    def genLine(self, event):
        if event.inaxes:
            global comps, selected_comps, P, allProps, opts
            x_click = event.xdata
            y_click = event.ydata
            x0n = np.array([x_click,y_click,1-x_click-y_click])
            x = RCM(comps,selected_comps,P,allProps,opts,x0n,2).x
            fig, ax = RCMplot(x, self.selected_comps, self.opts)
            
            if self.canvas:
                self.canvas.deleteLater()
                self.toolbar.deleteLater()
            
            self.canvas = FigureCanvas(fig)
            self.toolbar = CompactNavigationToolbar(self.canvas, self)
            self.coord_label = QLabel("", self.toolbar)
            self.toolbar.addWidget(self.coord_label)
            
            self.plot_layout = QVBoxLayout(self.plot_frame)
            self.plot_layout.addWidget(self.toolbar)
            self.plot_layout.addWidget(self.canvas)
            
            self.canvas.draw()

            self.canvas.mpl_connect('button_press_event', self.click_plot)

    def on_hover(self, event):
        if event.inaxes:
            self.coord_label.setText(f"x={event.xdata:.2f}, y={event.ydata:.2f}")
        else:
            self.coord_label.setText("")

    def auto_generate(self):
        self.plot_frame.hide()
        self.plot_frame = QFrame()
        self.main_layout.addWidget(self.plot_frame)
        self.canvas = None
        self.plot_figure()

    def clear_plot(self):
        self.plot_frame.hide()
        self.plot_frame = QFrame()
        self.main_layout.addWidget(self.plot_frame)
        self.canvas = None
        self.clear_figure()

    def click_plot(self, event):
        if event.inaxes:
            self.plot_frame.hide()
            self.plot_frame = QFrame()
            self.main_layout.addWidget(self.plot_frame)
            self.canvas = None
            self.genLine(event)

    def adjust_parameters(self):
        self.input_params_window = InputParamsWindow()
        self.input_params_window.show()

    def plot_options(self):
        self.plot_options_window = PlotOptsWindow()
        self.plot_options_window.show()

    def solver_options(self):
        self.solver_options_window = SolverOptsWindow()
        self.solver_options_window.show()

    def go_back(self):
        self.hide()
        self.new_simulation_window = NewSimulationWindow()
        self.new_simulation_window.show()

    def save_variables(self):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Variables", "", "RCM Files (*.rcm)")
        if file_path: 
            global P, comps, selected_comps, allProps, lmopts, opts, NRTL_aij, NRTL_bij, NRTL_cij, TcCel, Pc, omega, antoine_params, PLXANT_params
            data = {
                    'P': P,
                    'comps': comps,
                    'selected_comps': selected_comps,
                    'allProps': allProps,
                    'lmopts': lmopts,
                    'opts': opts,
                    'NRTL_aij': NRTL_aij,
                    'NRTL_bij': NRTL_bij,
                    'NRTL_cij': NRTL_cij,
                    'TcCel': TcCel,
                    'Pc': Pc,
                    'omega': omega,
                    'antoine_params': antoine_params,
                    'PLXANT_params': PLXANT_params
                    }
            with open(file_path, 'wb') as file:
                pickle.dump(data, file)

class InputParamsWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Input Parameters")
        self.create_widgets()

    def create_widgets(self):
        self.central_widget = QWidget()
        layout = QVBoxLayout(self.central_widget)

        self.antoine_params_button = QPushButton("Antoine")
        self.antoine_params_button.clicked.connect(self.antoine_params_btn)
        layout.addWidget(self.antoine_params_button)

        self.PLXANT_params_button = QPushButton("PLXANT")
        self.PLXANT_params_button.clicked.connect(self.PLXANT_params_btn)
        layout.addWidget(self.PLXANT_params_button)

        self.NRTL_aij_button = QPushButton("NRTL aij")
        self.NRTL_aij_button.clicked.connect(self.NRTL_aij_btn)
        layout.addWidget(self.NRTL_aij_button)

        self.NRTL_bij_button = QPushButton("NRTL bij")
        self.NRTL_bij_button.clicked.connect(self.NRTL_bij_btn)
        layout.addWidget(self.NRTL_bij_button)

        self.NRTL_cij_button = QPushButton("NRTL cij")
        self.NRTL_cij_button.clicked.connect(self.NRTL_cij_btn)
        layout.addWidget(self.NRTL_cij_button)

        self.SRK_params_button = QPushButton("SRK Parameters")
        self.SRK_params_button.clicked.connect(self.SRK_btn)
        layout.addWidget(self.SRK_params_button)
        
        self.setCentralWidget(self.central_widget)

    def antoine_params_btn(self):
        self.antoine_params_window = antoine_params_input()
        self.antoine_params_window.show() 

    def PLXANT_params_btn(self):
        self.PLXANT_params_window = PLXANT_params_input()
        self.PLXANT_params_window.show() 

    def NRTL_aij_btn(self):
        self.NRTL_aij_window = NRTL_aij_input()
        self.NRTL_aij_window.show() 

    def NRTL_bij_btn(self):
        self.NRTL_bij_window = NRTL_bij_input()
        self.NRTL_bij_window.show() 

    def NRTL_cij_btn(self):
        self.NRTL_cij_window = NRTL_cij_input()
        self.NRTL_cij_window.show() 

    def SRK_btn(self):
        self.SRK_window = SRK_input()
        self.SRK_window.show() 

class antoine_params_input(QWidget):
    def __init__(self):
        global comps, antoine_params 
        super().__init__()
        self.setGeometry(100, 100, 600, 400)
        self.antoine_bank = ['A', 'B', 'C']
        self.comps = comps
        self.antoine_params = antoine_params 
        self.create_widgets()
    
    def create_widgets(self):
        self.setWindowTitle('Antoine Parameters')

        self.table = CustomTableWidget(len(self.comps), len(self.antoine_bank))
        self.table.setHorizontalHeaderLabels([self.antoine_bank[i] for i in range(len(self.antoine_bank))])
        self.table.setVerticalHeaderLabels([self.comps[i] for i in range(len(self.comps))])

        self.populate_table()

        layout = QVBoxLayout()
        layout.addWidget(self.table)
        
        button_layout = QHBoxLayout()
        self.save_button = QPushButton('Save')
        self.save_button.clicked.connect(self.save_data)
        button_layout.addWidget(self.save_button)
        
        self.help_button = QPushButton('Help')
        self.help_button.clicked.connect(self.show_help)
        button_layout.addWidget(self.help_button)
        
        layout.addLayout(button_layout)
        self.setLayout(layout)

    def populate_table(self):
        for row in range(self.antoine_params.shape[0]):
            for col in range(self.antoine_params.shape[1]):
                self.table.setItem(row, col, QTableWidgetItem(str(self.antoine_params[row, col])))

    def save_data(self):
        global antoine_params 
        row_count = self.table.rowCount()
        col_count = self.table.columnCount()
        
        data = np.zeros((row_count, col_count))
        
        for row in range(row_count):
            for col in range(col_count):
                item = self.table.item(row, col)
                if item is not None and item.text() != '':
                    try:
                        data[row, col] = float(item.text())
                    except ValueError:
                        QMessageBox.warning(self, 'Invalid Input', f'Invalid value at row {row+1}, col {col+1}. Please enter numeric values.')
                        return
        
        antoine_params = data
    
    def show_help(self):
        QMessageBox.information(self, 'Help', 'Enter numeric values into the table. Click Save to save the data.')
        
class PLXANT_params_input(QWidget):
    def __init__(self):
        global comps, PLXANT_params 
        super().__init__()
        self.setGeometry(100, 100, 600, 400)
        self.PLXANT_bank = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']
        self.comps = comps
        self.PLXANT_params = PLXANT_params 
        self.create_widgets()
    
    def create_widgets(self):
        self.setWindowTitle('Extended Antoine Parameters')

        self.table = CustomTableWidget(len(self.comps), len(self.PLXANT_bank))
        self.table.setHorizontalHeaderLabels([self.PLXANT_bank[i] for i in range(len(self.PLXANT_bank))])
        self.table.setVerticalHeaderLabels([self.comps[i] for i in range(len(self.comps))])

        self.populate_table()

        layout = QVBoxLayout()
        layout.addWidget(self.table)
        
        button_layout = QHBoxLayout()
        self.save_button = QPushButton('Save')
        self.save_button.clicked.connect(self.save_data)
        button_layout.addWidget(self.save_button)
        
        self.help_button = QPushButton('Help')
        self.help_button.clicked.connect(self.show_help)
        button_layout.addWidget(self.help_button)
        
        layout.addLayout(button_layout)
        self.setLayout(layout)

    def populate_table(self):
        for row in range(self.PLXANT_params.shape[0]):
            for col in range(self.PLXANT_params.shape[1]):
                self.table.setItem(row, col, QTableWidgetItem(str(self.PLXANT_params[row, col])))

    def save_data(self):
        global PLXANT_params 
        row_count = self.table.rowCount()
        col_count = self.table.columnCount()
        
        data = np.zeros((row_count, col_count))
        
        for row in range(row_count):
            for col in range(col_count):
                item = self.table.item(row, col)
                if item is not None and item.text() != '':
                    try:
                        data[row, col] = float(item.text())
                    except ValueError:
                        QMessageBox.warning(self, 'Invalid Input', f'Invalid value at row {row+1}, col {col+1}. Please enter numeric values.')
                        return
        
        PLXANT_params = data
    
    def show_help(self):
        QMessageBox.information(self, 'Help', 'Enter numeric values into the table. Click Save to save the data.')
        
class NRTL_aij_input(QWidget):
    def __init__(self):
        global comps, NRTL_aij
        super().__init__()
        self.setGeometry(100, 100, 600, 400)
        self.comps = comps
        self.NRTL_aij = NRTL_aij
        self.create_widgets()
    
    def create_widgets(self):
        self.setWindowTitle('NRTL aij Parameters')

        self.table = CustomTableWidget(len(self.comps), len(self.comps))
        self.table.setHorizontalHeaderLabels([self.comps[i] for i in range(len(self.comps))])
        self.table.setVerticalHeaderLabels([self.comps[i] for i in range(len(self.comps))])

        self.populate_table()

        layout = QVBoxLayout()
        layout.addWidget(self.table)
        
        button_layout = QHBoxLayout()
        self.save_button = QPushButton('Save')
        self.save_button.clicked.connect(self.save_data)
        button_layout.addWidget(self.save_button)
        
        self.help_button = QPushButton('Help')
        self.help_button.clicked.connect(self.show_help)
        button_layout.addWidget(self.help_button)
        
        layout.addLayout(button_layout)
        self.setLayout(layout)

    def populate_table(self):
        for row in range(self.NRTL_aij.shape[0]):
            for col in range(self.NRTL_aij.shape[1]):
                self.table.setItem(row, col, QTableWidgetItem(str(self.NRTL_aij[row, col])))

    def save_data(self):
        global NRTL_aij
        row_count = self.table.rowCount()
        col_count = self.table.columnCount()
        
        data = np.zeros((row_count, col_count))
        
        for row in range(row_count):
            for col in range(col_count):
                item = self.table.item(row, col)
                if item is not None and item.text() != '':
                    try:
                        data[row, col] = float(item.text())
                    except ValueError:
                        QMessageBox.warning(self, 'Invalid Input', f'Invalid value at row {row+1}, col {col+1}. Please enter numeric values.')
                        return
        
        NRTL_aij = data
    
    def show_help(self):
        QMessageBox.information(self, 'Help', 'Enter numeric values into the table. Click Save to save the data.')

class NRTL_bij_input(QWidget):
    def __init__(self):
        global comps, NRTL_bij
        super().__init__()
        self.setGeometry(100, 100, 600, 400)
        self.comps = comps
        self.NRTL_bij = NRTL_bij
        self.create_widgets()
    
    def create_widgets(self):
        self.setWindowTitle('NRTL bij Parameters')

        self.table = CustomTableWidget(len(self.comps), len(self.comps))
        self.table.setHorizontalHeaderLabels([self.comps[i] for i in range(len(self.comps))])
        self.table.setVerticalHeaderLabels([self.comps[i] for i in range(len(self.comps))])

        self.populate_table()

        layout = QVBoxLayout()
        layout.addWidget(self.table)
        
        button_layout = QHBoxLayout()
        self.save_button = QPushButton('Save')
        self.save_button.clicked.connect(self.save_data)
        button_layout.addWidget(self.save_button)
        
        self.help_button = QPushButton('Help')
        self.help_button.clicked.connect(self.show_help)
        button_layout.addWidget(self.help_button)
        
        layout.addLayout(button_layout)
        self.setLayout(layout)

    def populate_table(self):
        for row in range(self.NRTL_bij.shape[0]):
            for col in range(self.NRTL_bij.shape[1]):
                self.table.setItem(row, col, QTableWidgetItem(str(self.NRTL_bij[row, col])))

    def save_data(self):
        global NRTL_bij
        row_count = self.table.rowCount()
        col_count = self.table.columnCount()
        
        data = np.zeros((row_count, col_count))
        
        for row in range(row_count):
            for col in range(col_count):
                item = self.table.item(row, col)
                if item is not None and item.text() != '':
                    try:
                        data[row, col] = float(item.text())
                    except ValueError:
                        QMessageBox.warning(self, 'Invalid Input', f'Invalid value at row {row+1}, col {col+1}. Please enter numeric values.')
                        return
        
        NRTL_bij = data
    
    def show_help(self):
        QMessageBox.information(self, 'Help', 'Enter numeric values into the table. Click Save to save the data.')

class NRTL_cij_input(QWidget):
    def __init__(self):
        global comps, NRTL_cij
        super().__init__()
        self.setGeometry(100, 100, 600, 400)
        self.comps = comps
        self.NRTL_cij = NRTL_cij
        self.create_widgets()
    
    def create_widgets(self):
        self.setWindowTitle('NRTL cij Parameters')

        self.table = CustomTableWidget(len(self.comps), len(self.comps))
        self.table.setHorizontalHeaderLabels([self.comps[i] for i in range(len(self.comps))])
        self.table.setVerticalHeaderLabels([self.comps[i] for i in range(len(self.comps))])

        self.populate_table()

        layout = QVBoxLayout()
        layout.addWidget(self.table)
        
        button_layout = QHBoxLayout()
        self.save_button = QPushButton('Save')
        self.save_button.clicked.connect(self.save_data)
        button_layout.addWidget(self.save_button)
        
        self.help_button = QPushButton('Help')
        self.help_button.clicked.connect(self.show_help)
        button_layout.addWidget(self.help_button)
        
        layout.addLayout(button_layout)
        self.setLayout(layout)

    def populate_table(self):
        for row in range(self.NRTL_cij.shape[0]):
            for col in range(self.NRTL_cij.shape[1]):
                self.table.setItem(row, col, QTableWidgetItem(str(self.NRTL_cij[row, col])))

    def save_data(self):
        global NRTL_cij
        row_count = self.table.rowCount()
        col_count = self.table.columnCount()
        
        data = np.zeros((row_count, col_count))
        
        for row in range(row_count):
            for col in range(col_count):
                item = self.table.item(row, col)
                if item is not None and item.text() != '':
                    try:
                        data[row, col] = float(item.text())
                    except ValueError:
                        QMessageBox.warning(self, 'Invalid Input', f'Invalid value at row {row+1}, col {col+1}. Please enter numeric values.')
                        return
        
        NRTL_cij = data
    
    def show_help(self):
        QMessageBox.information(self, 'Help', 'Enter numeric values into the table. Click Save to save the data.')

class SRK_input(QWidget):
    def __init__(self):
        global comps 
        super().__init__()
        self.setGeometry(100, 100, 600, 400)
        self.crit_params_bank = ['$T_C$', '$P_C$', '$\\omega$']
        self.comps = comps
        # self.crit_params = crit_params 
        self.create_widgets()
    
    def create_widgets(self):
        self.setWindowTitle('NRTL cij Parameters')

        self.table = CustomTableWidget(len(self.comps), len(self.crit_params_bank))
        self.table.setHorizontalHeaderLabels([self.crit_params_bank[i] for i in range(len(self.crit_params_bank))])
        self.table.setVerticalHeaderLabels([self.comps[i] for i in range(len(self.comps))])

        self.populate_table()

        layout = QVBoxLayout()
        layout.addWidget(self.table)
        
        button_layout = QHBoxLayout()
        self.save_button = QPushButton('Save')
        self.save_button.clicked.connect(self.save_data)
        button_layout.addWidget(self.save_button)
        
        self.help_button = QPushButton('Help')
        self.help_button.clicked.connect(self.show_help)
        button_layout.addWidget(self.help_button)
        
        layout.addLayout(button_layout)
        self.setLayout(layout)

    def populate_table(self):
        global TcCel, Pc, omega
        for row in range(len(TcCel)):
            self.table.setItem(row, 0, QTableWidgetItem(str(TcCel[row])))
        for row in range(len(Pc)):
            self.table.setItem(row, 1, QTableWidgetItem(str(Pc[row])))
        for row in range(len(omega)):
            self.table.setItem(row, 2, QTableWidgetItem(str(omega[row])))

    def save_data(self):
        global TcCel, Pc, omega 
        row_count = self.table.rowCount()
        col_count = self.table.columnCount()
        
        data = np.zeros((row_count, col_count))
        
        for row in range(row_count):
            for col in range(col_count):
                item = self.table.item(row, col)
                if item is not None and item.text() != '':
                    try:
                        data[row, col] = float(item.text())
                    except ValueError:
                        QMessageBox.warning(self, 'Invalid Input', f'Invalid value at row {row+1}, col {col+1}. Please enter numeric values.')
                        return
        
        TcCel = data[:,0]
        Pc = data[:,1]
        omega = data[:,2]
    
    def show_help(self):
        QMessageBox.information(self, 'Help', 'Enter numeric values into the table. Click Save to save the data.')

class PlotOptsWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Plot Options")
        self.create_widgets()

    def create_widgets(self):
        global opts

        self.central_widget = QWidget()
        layout = QGridLayout(self.central_widget)

        line_width_label = QLabel("Line Width:", self)
        layout.addWidget(line_width_label, 0, 0)

        self.linewidth_entry = QLineEdit(str(opts.linewidth), self)
        self.linewidth_entry.returnPressed.connect(self.save_linewidth)
        layout.addWidget(self.linewidth_entry, 0, 1)

        num_lines_label = QLabel("Auto-Gen # of Lines:", self)
        layout.addWidget(num_lines_label, 1, 0)

        self.num_lines_entry = QLineEdit(str(opts.lines), self)
        self.num_lines_entry.returnPressed.connect(self.save_num_lines)
        layout.addWidget(self.num_lines_entry, 1, 1)

        self.central_widget.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setCentralWidget(self.central_widget)

    def save_linewidth(self):
        global opts 
        try:
            opts.linewidth = float(self.linewidth_entry.text())
        except ValueError:
            QMessageBox.critical(self, "Invalid Input", "Please enter a valid number for line width.")

    def save_num_lines(self):
        global opts 
        try:
            opts.lines = int(self.num_lines_entry.text())
        except ValueError:
            QMessageBox.critical(self, "Invalid Input", "Please enter a valid interger for the amount of lines to genertate with auto-gen.")

class SolverOptsWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Solver Options")
        self.create_widgets()

    def create_widgets(self):
        global opts, lmopts

        self.central_widget = QWidget()
        layout = QGridLayout(self.central_widget)

        num_points_label = QLabel("Num. Points fwd/back:", self)
        layout.addWidget(num_points_label, 0, 0)

        self.num_points_entry = QLineEdit(str(opts.n_it), self)
        self.num_points_entry.returnPressed.connect(self.save_num_points)
        layout.addWidget(self.num_points_entry, 0, 1)

        maxiter_label = QLabel("Maximum Iterations:", self)
        layout.addWidget(maxiter_label, 1, 0)

        self.maxiter_entry = QLineEdit(str(lmopts['maxiter']), self)
        self.maxiter_entry.returnPressed.connect(self.save_maxiter)
        layout.addWidget(self.maxiter_entry, 1, 1)

        ftol_label = QLabel("Obj. Func. Tolerance:", self)
        layout.addWidget(ftol_label, 2, 0)

        self.ftol_entry = QLineEdit(str(lmopts['ftol']), self)
        self.ftol_entry.returnPressed.connect(self.save_ftol)
        layout.addWidget(self.ftol_entry, 2, 1)

        xtol_label = QLabel("Mol. Frac. Tolerance:", self)
        layout.addWidget(xtol_label, 3, 0)

        self.xtol_entry = QLineEdit(str(lmopts['xtol']), self)
        self.xtol_entry.returnPressed.connect(self.save_xtol)
        layout.addWidget(self.xtol_entry, 3, 1)
        
        self.central_widget.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setCentralWidget(self.central_widget)

    def save_num_points(self):
        global opts
        try:
            opts.n_it = int(self.num_points_entry.text())
        except ValueError:
            QMessageBox.critical(self, "Invalid Input", "Please enter a valid number for number of points to simulate forward and backwards.")

    def save_maxiter(self):
        global opts 
        try:
            lmopts['maxiter'] = int(self.maxiter_entry.text())
        except ValueError:
            QMessageBox.critical(self, "Invalid Input", "Please enter a valid integer for the amount of lines to genertate with auto-gen.")

    def save_ftol(self):
        global opts 
        try:
            lmopts['ftol'] = float(self.ftol_entry.text())
        except ValueError:
            QMessageBox.critical(self, "Invalid Input", "Please enter a valid number for the objective function tolerance.")

    def save_xtol(self):
        global opts 
        try:
            lmopts['xtol'] = float(self.xtol_entry.text())
        except ValueError:
            QMessageBox.critical(self, "Invalid Input", "Please enter a valid number for mole fraction tolerance.")

class CustomTableWidget(QTableWidget):
    def keyPressEvent(self, event):
        if event.modifiers() == Qt.KeyboardModifier.ControlModifier and event.key() == Qt.Key.Key_V:
            self.paste_data()
        else:
            super().keyPressEvent(event)
    
    def paste_data(self):
        clipboard = QApplication.clipboard()
        clipboard_text = clipboard.text()
        rows = clipboard_text.split('\n')
        
        for row_idx, row_data in enumerate(rows):
            columns = row_data.split('\t')
            for col_idx, col_data in enumerate(columns):
                if row_idx < self.rowCount() and col_idx < self.columnCount():
                    self.setItem(row_idx, col_idx, QTableWidgetItem(col_data))

class CompactNavigationToolbar(NavigationToolbar):
    toolitems = [t for t in NavigationToolbar.toolitems if
                 t[0] in ('Home', 'Pan', 'Zoom', 'Save')]

    def __init__(self, canvas, parent, coordinates=True):
        super().__init__(canvas, parent, coordinates)
        self.setIconSize(QSize(16, 16))  

class ClearFocusLineEdit(QLineEdit):
    def focusOutEvent(self, event):
        super().focusOutEvent(event)
        self.clearFocus()

if __name__ == "__main__":
    app = QApplication([])
    window = GetStartedWindow()
    window.show()
    app.exec()
