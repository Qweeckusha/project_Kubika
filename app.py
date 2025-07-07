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
        self.plot_resolution = 1000
        self.default_xlim = [-15, 15]
        self.default_ylim = [-15, 15]
        self.default_rlim = 15

        # <<< 1. Добавляем атрибут для хранения ID таймера >>>
        self.after_id = None

        self.create_widgets()
        self.set_strophoid_example()

    def create_widgets(self):
        # ... (код без изменений) ...
        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        ttk.Label(control_frame, text="Коэффициенты F(x,y)=0:", font=("Arial", 11, "bold")).grid(row=0, columnspan=2,
                                                                                                 pady=5)
        labels = ["A (x³)", "B (x²y)", "C (xy²)", "D (y³)", "E (x²)", "F (xy)", "G (y²)", "H (x)", "I (y)", "J"]
        for i, (key, var) in enumerate(self.coeffs.items()):
            ttk.Label(control_frame, text=f"{labels[i]}:").grid(row=i + 1, column=0, sticky='w', padx=5)
            ttk.Entry(control_frame, textvariable=var, width=8).grid(row=i + 1, column=1, pady=2)
        plot_button = ttk.Button(control_frame, text="Построить / Обновить", command=self.draw_plots_general)
        plot_button.grid(row=11, columnspan=2, pady=5, sticky='ew')
        reset_view_button = ttk.Button(control_frame, text="Сбросить вид", command=self.reset_views)
        reset_view_button.grid(row=12, columnspan=2, pady=(0, 5), sticky='ew')
        ttk.Separator(control_frame, orient='horizontal').grid(row=13, columnspan=2, sticky='ew', pady=10)
        ttk.Label(control_frame, text="Примеры:").grid(row=14, columnspan=2)
        strophoid_btn = ttk.Button(control_frame, text="Строфоида", command=self.set_strophoid_example)
        strophoid_btn.grid(row=15, columnspan=2, pady=2, sticky='ew')
        parabola_btn = ttk.Button(control_frame, text="Куб. парабола", command=self.set_parabola_example)
        parabola_btn.grid(row=16, columnspan=2, pady=2, sticky='ew')
        decartes_btn = ttk.Button(control_frame, text="Декартов лист", command=self.set_decartes_leaf_example)
        decartes_btn.grid(row=17, columnspan=2, pady=2, sticky='ew')
        plot_frame = ttk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.figure = plt.figure(figsize=(12, 6))
        self.ax1 = self.figure.add_subplot(1, 2, 1)
        self.ax2_polar = self.figure.add_subplot(1, 2, 2, projection='polar')
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)

    # <<< 2. Изменяем обработчик скролла >>>
    def on_scroll(self, event):
        ax = event.inaxes
        if ax is None: return

        # Сначала отменяем предыдущий запланированный вызов, если он был
        if self.after_id:
            self.root.after_cancel(self.after_id)

        # Немедленно обновляем только пределы осей, без перерисовки данных
        scale_factor = 1 / 1.2 if event.button == 'up' else 1.2

        if ax == self.ax1:
            cur_xlim = ax.get_xlim();
            cur_ylim = ax.get_ylim()
            xdata, ydata = event.xdata, event.ydata
            if xdata is None or ydata is None: return
            new_xlim = [(cur_xlim[0] - xdata) * scale_factor + xdata, (cur_xlim[1] - xdata) * scale_factor + xdata]
            new_ylim = [(cur_ylim[0] - ydata) * scale_factor + ydata, (cur_ylim[1] - ydata) * scale_factor + ydata]
            ax.set_xlim(new_xlim);
            ax.set_ylim(new_ylim)

            # Планируем "тяжелую" перерисовку через 150 мс
            self.after_id = self.root.after(150, self.update_cartesian_plot)

        elif ax == self.ax2_polar:
            cur_rlim = ax.get_ylim()
            new_r_max = cur_rlim[1] / scale_factor
            if new_r_max < 0.01: new_r_max = 0.01
            ax.set_ylim(0, new_r_max)

        # Быстрая перерисовка только осей и сетки, без данных
        self.canvas.draw_idle()

    # Остальные функции остаются без изменений
    def reset_views(self):
        self.ax1.set_xlim(self.default_xlim)
        self.ax1.set_ylim(self.default_ylim)
        self.update_cartesian_plot()
        self.ax2_polar.set_ylim(0, self.default_rlim)
        self.canvas.draw_idle()

    def draw_plots_general(self):
        self.draw_cartesian_contour()
        self.draw_polar_general()
        self.figure.tight_layout(pad=2.0)
        self.canvas.draw()

    def draw_cartesian_contour(self):
        self.ax1.clear()
        self.ax1.set_xlim(self.default_xlim)
        self.ax1.set_ylim(self.default_ylim)
        self.update_cartesian_plot()

    def update_cartesian_plot(self):
        self.after_id = None  # Сбрасываем ID, так как функция выполнилась
        xlim = self.ax1.get_xlim()
        ylim = self.ax1.get_ylim()

        # Очищаем только старый контур, а не все оси
        for collection in self.ax1.collections:
            collection.remove()

        c = {key: var.get() for key, var in self.coeffs.items()}

        x = np.linspace(xlim[0], xlim[1], self.plot_resolution)
        y = np.linspace(ylim[0], ylim[1], self.plot_resolution)
        X, Y = np.meshgrid(x, y)

        try:
            Z = (c['A'] * X ** 3 + c['B'] * X ** 2 * Y + c['C'] * X * Y ** 2 + c['D'] * Y ** 3 +
                 c['E'] * X ** 2 + c['F'] * X * Y + c['G'] * Y ** 2 + c['H'] * X + c['I'] * Y + c['J'])
        except MemoryError:
            print("Ошибка: не хватает памяти для построения графика с таким разрешением.")
            return

        self.ax1.contour(X, Y, Z, levels=[0], colors='blue', linewidths=1.5)
        self.ax1.set_title("Декартова система (общий метод)")
        self.ax1.set_xlabel("x")
        self.ax1.set_ylabel("y")
        self.ax1.grid(True)
        self.ax1.set_aspect('equal', adjustable='box')
        self.ax1.set_xlim(xlim)
        self.ax1.set_ylim(ylim)
        self.canvas.draw_idle()

    def draw_polar_general(self):
        ax = self.ax2_polar
        ax.clear()
        c = {key: var.get() for key, var in self.coeffs.items()}
        theta = np.linspace(0, 2 * np.pi, 4000)
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        P = c['A'] * cos_t ** 3 + c['B'] * cos_t ** 2 * sin_t + c['C'] * cos_t * sin_t ** 2 + c['D'] * sin_t ** 3
        Q = c['E'] * cos_t ** 2 + c['F'] * cos_t * sin_t + c['G'] * sin_t ** 2
        S = c['H'] * cos_t + c['I'] * sin_t
        T = c['J']
        r_branches = [[], [], []]
        for i in range(len(theta)):
            coeffs_poly = [P[i], Q[i], S[i], T]
            if np.all(np.isclose(coeffs_poly, 0)): continue
            roots = np.roots(coeffs_poly)
            real_roots = np.sort(roots[np.isreal(roots)].real)
            for j in range(3):
                r_branches[j].append(real_roots[j] if j < len(real_roots) else np.nan)
        ax.set_title("Полярная система (общий метод)")
        ax.grid(True)
        colors = ['purple', 'red', 'green']
        for i, r_branch in enumerate(r_branches):
            ax.plot(theta, r_branch, color=colors[i])
        ax.set_ylim(0, self.default_rlim)

    def clear_coeffs(self):
        for var in self.coeffs.values():
            var.set(0.0)

    def set_strophoid_example(self):
        self.clear_coeffs()
        a = 2.0
        self.coeffs['A'].set(1.0)
        self.coeffs['C'].set(1.0)
        self.coeffs['E'].set(-a)
        self.coeffs['G'].set(a)
        self.draw_plots_general()
        self.ax2_polar.set_ylim(0, a * 1.5)
        self.canvas.draw_idle()

    def set_parabola_example(self):
        self.clear_coeffs()
        self.coeffs['A'].set(0.2)
        self.coeffs['H'].set(-2.0)
        self.coeffs['I'].set(-1.0)
        self.coeffs['J'].set(1.0)
        self.draw_plots_general()
        self.ax2_polar.set_ylim(0, self.default_rlim)
        self.canvas.draw_idle()

    def set_decartes_leaf_example(self):
        self.clear_coeffs()
        a = 3.0
        self.coeffs['A'].set(1.0)
        self.coeffs['D'].set(1.0)
        self.coeffs['F'].set(-3.0 * a)
        self.draw_plots_general()
        self.ax2_polar.set_ylim(0, a * 1.5)
        self.canvas.draw_idle()


if __name__ == "__main__":
    root = tk.Tk()
    app = UniversalCubicPlotterApp(root)
    root.mainloop()