import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
import numpy as np


class UniversalCubicPlotterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Универсальный построитель кривых Кубика")
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        # --- Переменные состояния ---
        self.animation_var = tk.BooleanVar(value=False)
        self.animation = None
        self.anim_point_cartesian = None
        self.anim_point_polar = None
        self.path_data = None
        self.is_animation_running = False  # Флаг для блокировки зума

        # --- УПРАВЛЕНИЕ СКОРОСТЬЮ ---
        # Единый шаг для обеих систем координат.
        self.animation_frame_step = 14

        self.coeffs = {chr(65 + i): tk.DoubleVar(value=0.0) for i in range(10)}
        self.plot_resolution = 1000
        self.default_xlim = [-15, 15]
        self.default_ylim = [-15, 15]
        self.default_rlim = 15
        self.after_id = None

        self.create_widgets()
        self.set_strophoid_example()

    def on_closing(self):
        self.stop_animation()
        if self.after_id: self.root.after_cancel(self.after_id)
        plt.close('all')
        self.root.destroy()

    def create_widgets(self):
        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        ttk.Label(control_frame, text="Коэффициенты F(x,y)=0:", font=("Arial", 11, "bold")).grid(row=0, columnspan=2,
                                                                                                 pady=5)
        labels = ["A (x³)", "B (x²y)", "C (xy²)", "D (y³)", "E (x²)", "F (xy)", "G (y²)", "H (x)", "I (y)", "J"]
        for i, (key, var) in enumerate(self.coeffs.items()):
            ttk.Label(control_frame, text=f"{labels[i]}:").grid(row=i + 1, column=0, sticky='w', padx=5)
            ttk.Entry(control_frame, textvariable=var, width=8).grid(row=i + 1, column=1, pady=2)
        plot_button = ttk.Button(control_frame, text="Построить / Обновить", command=self.reset_views)
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
        ttk.Separator(control_frame, orient='horizontal').grid(row=18, columnspan=2, sticky='ew', pady=10)
        animation_check = ttk.Checkbutton(control_frame, text="Включить анимацию", variable=self.animation_var,
                                          command=self.toggle_animation)
        animation_check.grid(row=19, columnspan=2, sticky='w', pady=5)
        plot_frame = ttk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.figure = plt.figure(figsize=(12, 6))
        self.ax1 = self.figure.add_subplot(1, 2, 1)
        self.ax2_polar = self.figure.add_subplot(1, 2, 2, projection='polar')
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)

    def on_scroll(self, event):
        if self.is_animation_running: return
        ax = event.inaxes
        if ax is None: return
        if self.after_id:
            self.root.after_cancel(self.after_id)
            self.after_id = None

        if ax == self.ax1:
            scale_factor = 1.2 if event.button == 'down' else 1 / 1.2
            cur_xlim, cur_ylim = ax.get_xlim(), ax.get_ylim()
            xdata, ydata = event.xdata, event.ydata
            if xdata is None or ydata is None: return

            new_xlim = [(x - xdata) * scale_factor + xdata for x in cur_xlim]
            new_ylim = [(y - ydata) * scale_factor + ydata for y in cur_ylim]

            ax.set_xlim(new_xlim)
            ax.set_ylim(new_ylim)

            self.after_id = self.root.after(150, self.update_cartesian_plot)
            self.canvas.draw_idle()

        elif ax == self.ax2_polar:
            scale_factor = 1.2 if event.button == 'down' else 1 / 1.2
            cur_rlim = ax.get_ylim()
            new_r_max = cur_rlim[1] * scale_factor
            ax.set_ylim(0, max(0.01, new_r_max))
            self.canvas.draw_idle()

    def reset_views(self):
        was_running = self.is_animation_running
        self.stop_animation()
        self.draw_plots_general(reset_lims=True)
        if was_running:
            self.animation_var.set(True)
            self.start_animation()

    def clear_coeffs(self):
        for var in self.coeffs.values(): var.set(0.0)

    # --- Методы-примеры ---
    def set_strophoid_example(self):
        self.clear_coeffs();
        a = 2.0
        self.coeffs['A'].set(1.0);
        self.coeffs['C'].set(1.0);
        self.coeffs['E'].set(-a);
        self.coeffs['G'].set(a)
        self.draw_plots_general(reset_lims=True)

    def set_parabola_example(self):
        self.clear_coeffs()
        self.coeffs['A'].set(0.2);
        self.coeffs['H'].set(-2.0);
        self.coeffs['I'].set(-1.0);
        self.coeffs['J'].set(1.0)
        self.draw_plots_general(reset_lims=True)

    def set_decartes_leaf_example(self):
        self.clear_coeffs();
        a = 3.0
        self.coeffs['A'].set(1.0);
        self.coeffs['D'].set(1.0);
        self.coeffs['F'].set(-3.0 * a)
        self.draw_plots_general(reset_lims=True)

    # --- Основная логика отрисовки ---
    def draw_plots_general(self, reset_lims=False):
        self.stop_animation()
        self.update_cartesian_plot(reset_lims=reset_lims)
        self.draw_polar_general(reset_lims=reset_lims)
        self.figure.tight_layout(pad=2.0)
        self.canvas.draw()
        if self.animation_var.get():
            self.start_animation()

    def update_cartesian_plot(self, reset_lims=False):
        if self.after_id: self.root.after_cancel(self.after_id)

        xlim = self.default_xlim if reset_lims else self.ax1.get_xlim()
        ylim = self.default_ylim if reset_lims else self.ax1.get_ylim()

        for collection in self.ax1.collections:
            collection.remove()

        c = {key: var.get() for key, var in self.coeffs.items()}

        x_range = xlim[1] - xlim[0]
        y_range = ylim[1] - ylim[0]
        x_grid = np.linspace(xlim[0] - x_range * 0.1, xlim[1] + x_range * 0.1, self.plot_resolution)
        y_grid = np.linspace(ylim[0] - y_range * 0.1, ylim[1] + y_range * 0.1, self.plot_resolution)

        X, Y = np.meshgrid(x_grid, y_grid)
        try:
            Z = (c['A'] * X ** 3 + c['B'] * X ** 2 * Y + c['C'] * X * Y ** 2 + c['D'] * Y ** 3 +
                 c['E'] * X ** 2 + c['F'] * X * Y + c['G'] * Y ** 2 + c['H'] * X + c['I'] * Y + c['J'])
        except MemoryError:
            print("Ошибка: не хватает памяти.");
            return

        contour_set = self.ax1.contour(X, Y, Z, levels=[0], colors='blue', linewidths=1.5)

        if reset_lims:
            self.path_data = None
            if contour_set.allsegs and contour_set.allsegs[0]:
                all_vertices = np.vstack(contour_set.allsegs[0])
                x, y = all_vertices[:, 0], all_vertices[:, 1]
                self.path_data = {
                    'cartesian': all_vertices,
                    'polar': np.column_stack((np.arctan2(y, x), np.sqrt(x ** 2 + y ** 2)))
                }

        self.ax1.set_title("Декартова система (общий метод)")
        self.ax1.grid(True)
        self.ax1.set_aspect('equal', adjustable='box')
        self.ax1.set_xlim(xlim)
        self.ax1.set_ylim(ylim)

    def draw_polar_general(self, reset_lims=False):
        self.ax2_polar.clear()
        c = {key: var.get() for key, var in self.coeffs.items()}
        theta = np.linspace(0, 2 * np.pi, 4000)
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        P = c['A'] * cos_t ** 3 + c['B'] * cos_t ** 2 * sin_t + c['C'] * cos_t * sin_t ** 2 + c['D'] * sin_t ** 3
        Q = c['E'] * cos_t ** 2 + c['F'] * cos_t * sin_t + c['G'] * sin_t ** 2
        S = c['H'] * cos_t + c['I'] * sin_t
        T = c['J']

        r_branches = [[] for _ in range(3)]
        all_positive_r_values = []

        for i in range(len(theta)):
            coeffs = [P[i], Q[i], S[i], T]
            if np.all(np.isclose(coeffs, 0)): continue
            roots = np.roots(coeffs)
            real_roots = np.sort(roots[np.isreal(roots)].real)

            all_positive_r_values.extend(r for r in real_roots if (r >= 0 and np.isfinite(r) and r < 1000))

            for j in range(3):
                value_to_add = real_roots[j] if (j < len(real_roots) and real_roots[j] >= 0) else np.nan
                r_branches[j].append(value_to_add)

        colors = ['purple', 'red', 'green']
        for i, r_branch in enumerate(r_branches):
            self.ax2_polar.plot(theta, np.array(r_branch), color=colors[i], linewidth=1.5)

        self.ax2_polar.set_title("Полярная система (общий метод)")
        self.ax2_polar.grid(True)

        if reset_lims:
            if all_positive_r_values:
                max_r = max(all_positive_r_values)
                final_r_lim = np.clip(max_r * 1.2, 5, 25)
                self.ax2_polar.set_ylim(0, final_r_lim)
            else:
                self.ax2_polar.set_ylim(0, self.default_rlim)

    # --- Управление анимацией ---
    def toggle_animation(self):
        if self.animation_var.get():
            if not self.animation:
                self.start_animation()
            else:
                self.animation.resume(); self.is_animation_running = True
        elif self.animation:
            self.animation.pause(); self.is_animation_running = False; self.draw_plots_general(reset_lims=True)

    def start_animation(self):
        self.stop_animation()
        if not self.path_data: print("Нет данных для построения анимации."); return

        self.anim_point_cartesian, = self.ax1.plot([], [], 'ro', markersize=6)
        self.anim_point_polar, = self.ax2_polar.plot([], [], 'ro', markersize=6)

        self.animation = FuncAnimation(self.figure, self._animate_frame,
                                       init_func=self._init_animation, frames=len(self.path_data['cartesian']),
                                       interval=10, blit=True, repeat=True)

        self.is_animation_running = True
        self.canvas.draw()

    def stop_animation(self):

        if self.animation:
            self.animation.event_source.stop()
            self.animation = None

        for point_attr in ['anim_point_cartesian', 'anim_point_polar']:
            point = getattr(self, point_attr, None)
            if point:
                point.set_visible(False)
                setattr(self, point_attr, None)
        self.is_animation_running = False
        self.canvas.draw_idle()

    def _init_animation(self):
        self.anim_point_cartesian.set_data([], [])
        self.anim_point_polar.set_data([], [])
        return self.anim_point_cartesian, self.anim_point_polar

    def _animate_frame(self, i):
        path_index = (i * self.animation_frame_step) % len(self.path_data['cartesian'])

        cart_point = self.path_data['cartesian'][path_index]
        polar_point = self.path_data['polar'][path_index]

        self.anim_point_cartesian.set_data([cart_point[0]], [cart_point[1]])
        self.anim_point_polar.set_data([polar_point[0]], [polar_point[1]])

        return self.anim_point_cartesian, self.anim_point_polar


if __name__ == "__main__":
    root = tk.Tk()
    app = UniversalCubicPlotterApp(root)
    root.mainloop()