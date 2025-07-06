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
        self.current_curve_type = "cubic_parabola"

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

        control_frame = ttk.Frame(self.root)
        control_frame.pack(padx=10, pady=10, fill=tk.X)

        param_frame = ttk.Frame(control_frame)
        param_frame.pack(fill=tk.X, expand=True)
        ttk.Label(param_frame, text="Параметр a:").pack(side=tk.LEFT, padx=(0, 5))
        entry = ttk.Entry(param_frame, textvariable=self.a_var, width=8)
        entry.pack(side=tk.LEFT)
        entry.bind("<Return>", lambda event: self.redraw_current_curve())

        slider = ttk.Scale(param_frame, from_=-100.0, to=100.0, orient=tk.HORIZONTAL, variable=self.a_var,
                           command=self.on_slider_move)
        slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=10)

        action_frame = ttk.Frame(control_frame)
        action_frame.pack(fill=tk.X, pady=(5, 0))
        self.build_button = ttk.Button(action_frame, text="Построить", command=self.redraw_current_curve)
        self.build_button.pack(side=tk.LEFT)

        zoom_frame = ttk.Frame(action_frame)
        zoom_frame.pack(side=tk.LEFT, padx=30)
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

    def on_slider_move(self, value):
        self.a_var.set(round(float(value), 2))
        self.redraw_current_curve()

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
        self.current_curve_type = "cubic_parabola"; self.draw_cubic_parabola()

    def select_strophoid(self):
        self.current_curve_type = "strophoid"; self.draw_strophoid()

    def select_decartes_leaf(self):
        self.current_curve_type = "decartes_leaf"; self.draw_decartes_leaf()

    def redraw_current_curve(self):
        draw_func_name = f"draw_{self.current_curve_type}"
        draw_func = getattr(self, draw_func_name)
        draw_func()

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
        self._prepare_axes("Кубическая парабола", "Полярный аналог")
        self.ax1.set_aspect('auto')

        if a > 0:
            x_limit = np.sqrt(a) * 1.5 if a > 1 else 1.5
        else:
            x_limit = 1.5
        x = np.linspace(-x_limit, x_limit, 400)
        y = x ** 3 - a * x
        sign = "-" if a >= 0 else "+"
        self.ax1.plot(x, y, label=f"$y = x^3 {sign} {abs(a):.1f}x$")
        self.ax1.legend()
        self.ax1.set_xlim(-x_limit, x_limit)

        theta = np.linspace(0, 4 * np.pi, 2000)
        r = theta ** 3 - a * theta
        self.ax2.plot(theta, r)

        self.canvas.draw()

    def draw_strophoid(self):
        a = self.a_var.get()
        self._prepare_axes("Строфоида", "Строфоида (полярн.)")
        self.ax1.set_aspect('equal', 'box')
        t = np.linspace(-30, 30, 4000)
        x = a * (t ** 2 - 1) / (t ** 2 + 1)
        y = a * t * (t ** 2 - 1) / (t ** 2 + 1)
        self.ax1.plot(x, y)
        self.ax1.set_xlim(-abs(a) * 1.5, abs(a) * 1.5)
        self.ax1.set_ylim(-abs(a) * 1.5, abs(a) * 1.5)

        theta = np.linspace(0, 2 * np.pi, 2000)
        r = -a * np.cos(2 * theta) / np.cos(theta)
        r[np.abs(r) > 4 * abs(a)] = np.nan
        self.ax2.plot(theta, r)
        self.ax2.set_ylim(0, abs(a) * 2.5)
        self.canvas.draw()

    def draw_decartes_leaf(self):
        a = self.a_var.get()
        self._prepare_axes("Декартов лист", "Декартов лист (полярн.)")
        self.ax1.set_aspect('equal', 'box')

        theta = np.linspace(0, 2 * np.pi, 4000)

        with np.errstate(divide='ignore', invalid='ignore'):
            sin_t, cos_t = np.sin(theta), np.cos(theta)
            denominator = sin_t ** 3 + cos_t ** 3
            r = (3 * a * sin_t * cos_t) / denominator

        # Генерируем декартовы точки из полярных
        x = r * np.cos(theta)
        y = r * np.sin(theta)

        self.ax1.plot(x, y, color='blue')
        if a != 0:
            # Определяем границы для линии асимптоты
            xlim = self.ax1.get_xlim()
            asymptote_x = np.array(xlim)
            asymptote_y = -asymptote_x - a
            # Рисуем красную пунктирную линию
            self.ax1.plot(asymptote_x, asymptote_y, color='red', label='Асимптота')
            self.ax1.legend()
        self.ax1.set_xlim(-abs(a) * 1.5, abs(a) * 1.5)
        self.ax1.set_ylim(-abs(a) * 1.5, abs(a) * 1.5)

        # Полярный график строим по той же самой, правильной формуле
        self.ax2.plot(theta, r)
        self.ax2.set_ylim(0, abs(a) * 1.5)

        self.canvas.draw()

    def on_resize(self, event):
        self.figure.tight_layout()
        if hasattr(self, 'canvas'): self.canvas.draw()


if __name__ == "__main__":
    root = tk.Tk()
    app = CubicPlotterApp(root)
    root.mainloop()