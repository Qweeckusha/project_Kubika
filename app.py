import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np


class UniversalCubicPlotterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Универсальный построитель кривых Кубика")

        self.coeffs = {chr(65 + i): tk.DoubleVar(value=0.0) for i in range(10)}
        self.create_widgets()
        self.set_strophoid_example()

    def create_widgets(self):
        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)

        ttk.Label(control_frame, text="Коэффициенты F(x,y)=0:", font=("Arial", 11, "bold")).grid(row=0, columnspan=2,
                                                                                                 pady=5)
        labels = ["A (x³)", "B (x²y)", "C (xy²)", "D (y³)", "E (x²)", "F (xy)", "G (y²)", "H (x)", "I (y)", "J"]
        for i, (key, var) in enumerate(self.coeffs.items()):
            ttk.Label(control_frame, text=f"{labels[i]}:").grid(row=i + 1, column=0, sticky='w', padx=5)
            ttk.Entry(control_frame, textvariable=var, width=8).grid(row=i + 1, column=1, pady=2)

        plot_button = ttk.Button(control_frame, text="Построить (общий метод)", command=self.draw_plots_general)
        plot_button.grid(row=11, columnspan=2, pady=10)

        ttk.Separator(control_frame, orient='horizontal').grid(row=12, columnspan=2, sticky='ew', pady=10)
        ttk.Label(control_frame, text="Примеры (точный метод):").grid(row=13, columnspan=2)
        strophoid_btn = ttk.Button(control_frame, text="Строфоида", command=self.set_strophoid_example)
        strophoid_btn.grid(row=14, columnspan=2, pady=2, sticky='ew')
        parabola_btn = ttk.Button(control_frame, text="Куб. парабола", command=self.set_parabola_example)
        parabola_btn.grid(row=15, columnspan=2, pady=2, sticky='ew')
        decartes_btn = ttk.Button(control_frame, text="Декартов лист", command=self.set_decartes_leaf_example)
        decartes_btn.grid(row=16, columnspan=2, pady=2, sticky='ew')

        plot_frame = ttk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.figure = plt.figure(figsize=(12, 6))
        self.ax1 = self.figure.add_subplot(1, 2, 1)
        self.ax2_polar = self.figure.add_subplot(1, 2, 2, projection='polar')

        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)

    # --- Общий метод отрисовки ---
    def draw_plots_general(self):
        self.draw_cartesian_contour()
        self.draw_polar_general()
        self.figure.tight_layout(pad=2.0)
        self.canvas.draw()

    def draw_cartesian_contour(self):
        # ... (код для contour, без изменений)
        ax = self.ax1
        ax.clear()
        c = {key: var.get() for key, var in self.coeffs.items()}

        x = np.linspace(-500, 500, 1000)
        y = np.linspace(-500, 500, 1000)
        X, Y = np.meshgrid(x, y)

        Z = (c['A'] * X ** 3 + c['B'] * X ** 2 * Y + c['C'] * X * Y ** 2 + c['D'] * Y ** 3 +
             c['E'] * X ** 2 + c['F'] * X * Y + c['G'] * Y ** 2 + c['H'] * X + c['I'] * Y + c['J'])

        ax.contour(X, Y, Z, levels=[0], colors='blue', linewidths=2)

        ax.set_title("Декартова система (общий метод)")
        ax.set_xlabel("x");
        ax.set_ylabel("y")
        ax.grid(True);
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim([-15, 15]);
        ax.set_ylim([-15, 15])

    def draw_polar_general(self):
        # ... (код для общего полярного графика, без изменений)
        ax = self.ax2_polar
        ax.clear()
        c = {key: var.get() for key, var in self.coeffs.items()}
        theta = np.linspace(0, 2 * np.pi, 2000)

        cos_t, sin_t = np.cos(theta), np.sin(theta)
        P = c['A'] * cos_t ** 3 + c['B'] * cos_t ** 2 * sin_t + c['C'] * cos_t * sin_t ** 2 + c['D'] * sin_t ** 3
        Q = c['E'] * cos_t ** 2 + c['F'] * cos_t * sin_t + c['G'] * sin_t ** 2
        S = c['H'] * cos_t + c['I'] * sin_t
        T = c['J']

        r_branches = [[], [], []]
        for i in range(len(theta)):
            coeffs = [P[i], Q[i], S[i], T]
            if not np.any(np.nan_to_num(coeffs)): continue
            roots = np.roots(coeffs)
            real_roots = np.sort(roots[np.isreal(roots)].real)
            for j in range(3):
                r_branches[j].append(real_roots[j] if j < len(real_roots) else np.nan)

        ax.set_title("Полярная система (общий метод)")
        ax.grid(True)
        colors = ['green', 'red', 'purple']
        for i, r_branch in enumerate(r_branches):
            ax.plot(theta, r_branch, color=colors[i])
        ax.set_ylim(bottom=0, top=15)

    # --- Специализированные методы отрисовки ---
    def draw_strophoid_parametric(self, a=2.0):
        # ... (код для строфоиды, без изменений)
        self.ax1.clear();
        self.ax2_polar.clear()
        t = np.linspace(-30, 30, 4000)
        x = a * (t ** 2 - 1) / (t ** 2 + 1)
        y = a * t * (t ** 2 - 1) / (t ** 2 + 1)
        self.ax1.plot(x, y, color='blue')
        self.ax1.set_title("Строфоида (точный метод)")
        self.ax1.grid(True);
        self.ax1.set_aspect('equal', 'box')
        self.ax1.set_xlim([-a * 1.5, a * 1.5]);
        self.ax1.set_ylim([-a * 1.5, a * 1.5])

        theta = np.linspace(0, 2 * np.pi, 2000)
        r = -a * np.cos(2 * theta) / np.cos(theta)
        r[np.abs(r) > 4 * a] = np.nan
        self.ax2_polar.plot(theta, r, color='green')
        self.ax2_polar.set_title("Строфоида (точный метод)")
        self.ax2_polar.grid(True);
        self.ax2_polar.set_ylim(bottom=0, top=2.5 * a)

        self.figure.tight_layout(pad=2.0)
        self.canvas.draw()

    def draw_parabola_explicit(self, coeffs):
        self.ax1.clear();
        self.ax2_polar.clear()

        # Декартов график
        x = np.linspace(-10, 10, 1000)
        y = coeffs['A'] * x ** 3 + coeffs['E'] * x ** 2 + coeffs['H'] * x + coeffs['J']
        self.ax1.plot(x, y, color='blue')
        self.ax1.set_title("Кубическая парабола (точный метод)")
        self.ax1.grid(True);
        self.ax1.set_ylim([-10, 10])

        # Полярный график
        theta = np.linspace(0, 2 * np.pi, 1000)
        r = coeffs['A'] * theta ** 3 + coeffs['E'] * theta ** 2 + coeffs['H'] * theta + coeffs['J']
        self.ax2_polar.plot(theta, r, color='green')
        self.ax2_polar.set_title("Кубическая парабола (точный метод)")
        self.ax2_polar.grid(True);
        self.ax2_polar.set_ylim([-10, 10])

        self.figure.tight_layout(pad=2.0)
        self.canvas.draw()

    def draw_decartes_leaf_parametric(self, a=3.0):
        self.ax1.clear();
        self.ax2_polar.clear()

        # Декартов график
        t = np.linspace(-50, 50, 4000)
        # Исключаем точку t=-1, где знаменатель равен нулю
        t = t[t != -1]
        x = (3 * a * t) / (1 + t ** 3)
        y = (3 * a * t ** 2) / (1 + t ** 3)

        self.ax1.plot(x, y, color='blue')
        self.ax1.set_title("Декартов лист (точный метод)")
        self.ax1.grid(True);
        self.ax1.set_aspect('equal', 'box')
        self.ax1.set_xlim([-2 * a, 2 * a]);
        self.ax1.set_ylim([-2 * a, 2 * a])

        # Полярный график
        theta = np.linspace(0, 2 * np.pi, 2000)
        denominator = np.sin(theta) ** 3 + np.cos(theta) ** 3
        # Исключаем деление на ноль
        r = np.full_like(theta, np.nan)
        mask = np.abs(denominator) > 1e-9
        r[mask] = (3 * a * np.sin(theta[mask]) * np.cos(theta[mask])) / denominator[mask]

        self.ax2_polar.plot(theta, r, color='green')
        self.ax2_polar.set_title("Декартов лист (точный метод)")
        self.ax2_polar.grid(True);
        self.ax2_polar.set_ylim(bottom=0, top=a * 1.5)

        self.figure.tight_layout(pad=2.0)
        self.canvas.draw()

    def on_scroll(self, event):
        # ... (код для скролла, без изменений)
        if event.inaxes != self.ax1: return
        ax = self.ax1
        scale_factor = 1.2 if event.button == 'up' else 1 / 1.2
        cur_xlim = ax.get_xlim();
        cur_ylim = ax.get_ylim()
        xdata, ydata = event.xdata, event.ydata
        new_xlim = [(cur_xlim[0] - xdata) * scale_factor + xdata, (cur_xlim[1] - xdata) * scale_factor + xdata]
        new_ylim = [(cur_ylim[0] - ydata) * scale_factor + ydata, (cur_ylim[1] - ydata) * scale_factor + ydata]
        ax.set_xlim(new_xlim);
        ax.set_ylim(new_ylim)
        self.canvas.draw_idle()

    # --- Функции-обработчики для кнопок ---
    def clear_coeffs(self):
        for var in self.coeffs.values(): var.set(0.0)

    def set_strophoid_example(self):
        self.clear_coeffs()
        a = 2.0
        self.coeffs['A'].set(1.0);
        self.coeffs['C'].set(1.0)
        self.coeffs['E'].set(-a);
        self.coeffs['G'].set(a)
        self.draw_strophoid_parametric(a)

    def set_parabola_example(self):
        self.clear_coeffs()
        # Для y = 0.2x³ - 2x + 1
        # Нужно преобразовать в y - (0.2x³ - 2x + 1) = 0
        # -0.2x³ + 2x + y - 1 = 0
        coeffs = {'A': 0.2, 'E': 0, 'H': -2.0, 'J': 1.0}
        self.coeffs['A'].set(coeffs['A'])
        self.coeffs['H'].set(coeffs['H'])
        self.coeffs['I'].set(-1.0)  # y -> -1*y
        self.coeffs['J'].set(coeffs['J'])
        self.draw_parabola_explicit(coeffs)

    def set_decartes_leaf_example(self):
        self.clear_coeffs()
        a = 3.0
        # x³ + y³ - 3axy = 0
        self.coeffs['A'].set(1.0)
        self.coeffs['D'].set(1.0)
        self.coeffs['F'].set(-3.0 * a)
        self.draw_decartes_leaf_parametric(a)


if __name__ == "__main__":
    root = tk.Tk()
    app = UniversalCubicPlotterApp(root)
    root.mainloop()