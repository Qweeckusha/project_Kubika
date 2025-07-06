import tkinter as tk
from tkinter import ttk
from tkinter import Menu
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np


class CubicPlotterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Визуализатор кривых Кубика")

        self.a_var = tk.DoubleVar(value=2.0)
        self.current_curve_type = "parabola"

        self.create_widgets()
        self.select_cubic_parabola()

    def create_widgets(self):
        menubar = Menu(self.root)
        self.root.config(menu=menubar)

        curve_menu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Выбор кривой", menu=curve_menu)

        curve_menu.add_command(label="Кубическая парабола", command=self.select_cubic_parabola)
        curve_menu.add_command(label="Строфоида", command=self.select_strophoid)
        curve_menu.add_command(label="Декартов лист", command=self.select_decartes_leaf)
        curve_menu.add_separator()
        curve_menu.add_command(label="Выход", command=self.root.quit)

        frame = ttk.Frame(self.root)
        frame.pack(padx=10, pady=10, fill=tk.X)

        ttk.Label(frame, text="Параметр a:").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(frame, textvariable=self.a_var, width=10).pack(side=tk.LEFT)

        self.build_button = ttk.Button(frame, text="Построить", command=self.redraw_current_curve)
        self.build_button.pack(side=tk.LEFT, padx=10)

        zoom_frame = ttk.Frame(frame)
        zoom_frame.pack(side=tk.LEFT, padx=20)
        ttk.Label(zoom_frame, text="Масштаб:").pack(side=tk.LEFT)
        zoom_in_button = ttk.Button(zoom_frame, text="+", width=3, command=lambda: self.zoom(0.8))
        zoom_in_button.pack(side=tk.LEFT)
        zoom_out_button = ttk.Button(zoom_frame, text="-", width=3, command=lambda: self.zoom(1.25))
        zoom_out_button.pack(side=tk.LEFT)
        reset_zoom_button = ttk.Button(zoom_frame, text="Сброс", command=self.redraw_current_curve)
        reset_zoom_button.pack(side=tk.LEFT, padx=5)

        self.figure = plt.figure(figsize=(12, 6))
        gs = GridSpec(1, 2, width_ratios=[3, 2], figure=self.figure)
        self.ax1 = self.figure.add_subplot(gs[0, 0])
        self.ax2 = self.figure.add_subplot(gs[0, 1], projection='polar')

        self.canvas = FigureCanvasTkAgg(self.figure, master=self.root)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.root.bind("<Configure>", self.on_resize)

    def _prepare_axes(self, cartesian_title, polar_title):
        self.ax1.clear()
        self.ax2.clear()
        self.ax1.set_title(cartesian_title)
        self.ax2.set_title(polar_title)
        self.ax1.grid(True)
        self.ax2.grid(True)
        self.ax1.set_xlabel("x")
        self.ax1.set_ylabel("y")

    def select_cubic_parabola(self):
        self.current_curve_type = "parabola"
        self.draw_cubic_parabola()

    def select_strophoid(self):
        self.current_curve_type = "strophoid"
        self.draw_strophoid()

    def select_decartes_leaf(self):
        self.current_curve_type = "leaf"
        self.draw_decartes_leaf()

    def redraw_current_curve(self):
        if self.current_curve_type == "parabola":
            self.draw_cubic_parabola()
        elif self.current_curve_type == "strophoid":
            self.draw_strophoid()
        elif self.current_curve_type == "leaf":
            self.draw_decartes_leaf()

    def zoom(self, factor):
        cur_xlim = self.ax1.get_xlim()
        cur_ylim = self.ax1.get_ylim()
        x_center = (cur_xlim[1] + cur_xlim[0]) / 2
        y_center = (cur_ylim[1] + cur_ylim[0]) / 2
        new_width = (cur_xlim[1] - cur_xlim[0]) * factor
        new_height = (cur_ylim[1] - cur_ylim[0]) * factor
        self.ax1.set_xlim([x_center - new_width / 2, x_center + new_width / 2])
        self.ax1.set_ylim([y_center - new_height / 2, y_center + new_height / 2])
        self.canvas.draw_idle()

    def draw_cubic_parabola(self):
        a = self.a_var.get()
        if a <= 0: a = 0.01
        self._prepare_axes("Кубическая парабола", "Полярный аналог")

        # <<< ГЛАВНЫЙ ФИКС: Явно сбрасываем соотношение сторон к автоматическому >>>
        self.ax1.set_aspect('auto')

        if a > 0:
            x_extremum = np.sqrt(a / 3)
            y_extremum = abs(x_extremum ** 3 - a * x_extremum)
            x_limit = x_extremum * 2.5
            y_limit = y_extremum * 1.5 if y_extremum > 0 else 10
        else:
            x_limit, y_limit = 4, 15

        x = np.linspace(-x_limit, x_limit, 400)
        y = x ** 3 - a * x
        self.ax1.plot(x, y, label=f"$y = x^3 - {a:.1f}x$")
        self.ax1.legend()
        self.ax1.set_xlim(-x_limit, x_limit)
        self.ax1.set_ylim(-y_limit, y_limit)

        theta = np.linspace(0, 2 * np.pi, 400)
        r = theta ** 3 - a * theta
        self.ax2.plot(theta, r)
        self.canvas.draw()

    def draw_strophoid(self):
        a = self.a_var.get()
        self._prepare_axes("Строфоида", "Строфоида (полярн.)")

        self.ax1.set_aspect('equal', 'box')  # Устанавливаем для этой кривой
        t = np.linspace(-30, 30, 4000)
        x = a * (t ** 2 - 1) / (t ** 2 + 1)
        y = a * t * (t ** 2 - 1) / (t ** 2 + 1)
        self.ax1.plot(x, y)
        self.ax1.set_xlim(-a * 1.5, a * 1.5);
        self.ax1.set_ylim(-a * 1.5, a * 1.5)

        theta = np.linspace(0, 2 * np.pi, 2000)
        r = -a * np.cos(2 * theta) / np.cos(theta)
        r[np.abs(r) > 4 * a] = np.nan
        self.ax2.plot(theta, r)
        self.ax2.set_ylim(0, 2.5 * a)
        self.canvas.draw()

    def draw_decartes_leaf(self):
        a = self.a_var.get()
        self._prepare_axes("Декартов лист", "Декартов лист (полярн.)")

        self.ax1.set_aspect('equal', 'box')  # Устанавливаем для этой кривой
        theta = np.linspace(0, 2 * np.pi, 4000)

        with np.errstate(divide='ignore', invalid='ignore'):
            sin_t, cos_t = np.sin(theta), np.cos(theta)
            denominator = sin_t ** 3 + cos_t ** 3
            r = (3 * a * sin_t * cos_t) / denominator

        x = r * np.cos(theta)
        y = r * np.sin(theta)

        self.ax1.plot(x, y, color='blue')
        self.ax1.set_xlim(-a * 1.5, a * 1.5);
        self.ax1.set_ylim(-a * 1.5, a * 1.5)

        self.ax2.plot(theta, r)
        self.ax2.set_ylim(0, a * 1.5)

        self.canvas.draw()

    def on_resize(self, event):
        self.figure.tight_layout()
        if hasattr(self, 'canvas'):
            self.canvas.draw()


if __name__ == "__main__":
    root = tk.Tk()
    app = CubicPlotterApp(root)
    root.mainloop()