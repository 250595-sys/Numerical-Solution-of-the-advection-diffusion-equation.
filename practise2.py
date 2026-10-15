import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button

# ==========================================
# 1. PARAMETERS & INITIAL CONDITION
# ==========================================
v = 1.0
del_t = 0.02
del_x = 0.05
alpha = v * del_t / del_x  # Courant number = 0.4
Nt = 100
Nx = 100

x = np.linspace(0, (Nx - 1) * del_x, Nx)
u0 = np.exp(-((x - 1.5) ** 2) / 0.5)  # Initial pulse at t=0


# ==========================================
# 2. NUMERICAL SCHEME FUNCTIONS
# ==========================================


# 1. FTCS (Forward Time Centered Space)
def solve_ftcs(u0, alpha, Nt, Nx):
  U = np.zeros((Nt, Nx))
  U[0, :] = u0.copy()
  for i in range(0, Nt - 1):
    for j in range(1, Nx - 1):
      U[i + 1, j] = U[i, j] - 0.5 * alpha * (U[i, j + 1] - U[i, j - 1])
      # Outflow (Transmissive) Boundary Condition
    U[i + 1, 0] = 0.0 # Inflow (left)
    U[i + 1, -1] = U[i + 1, -2] # Outflow (right)
  return U


# 2. Upwind Method
def solve_upwind(u0, alpha, Nt, Nx):
  U = np.zeros((Nt, Nx))
  U[0, :] = u0.copy()
  for i in range(0, Nt - 1):
    for j in range(1, Nx - 1):
      U[i + 1, j] = U[i, j] - alpha * (U[i, j] - U[i, j - 1])
    U[i + 1, 0] = 0.0
    U[i + 1, -1] = U[i + 1, -2]
  return U


# 3. Lax-Friedrichs
def solve_lax_friedrichs(u0, alpha, Nt, Nx):
  U = np.zeros((Nt, Nx))
  U[0, :] = u0.copy()
  for i in range(0, Nt - 1):
    for j in range(1, Nx - 1):
      U[i + 1, j] = 0.5 * (U[i, j + 1] + U[i, j - 1]) - 0.5 * alpha * (U[i, j + 1] - U[i, j - 1])
    U[i + 1, 0] = 0.0
    U[i + 1, -1] = U[i + 1, -2]
  return U


# 4. Lax-Wendroff
def solve_lax_wendroff(u0, alpha, Nt, Nx):
  U = np.zeros((Nt, Nx))
  U[0, :] = u0.copy()
  for i in range(0, Nt - 1):
    for j in range(1, Nx - 1):
      U[i + 1, j] = (U[i, j] - 0.5 * alpha * (U[i, j + 1] - U[i, j - 1]) + 0.5 * (alpha**2) * (U[i, j + 1] - 2 * U[i, j] + U[i, j - 1]))
    U[i + 1, 0] = 0.0
    U[i + 1, -1] = U[i + 1, -2]
  return U


# 5. Leapfrog (CTCS)
def solve_leapfrog(u0, alpha, Nt, Nx):
  U = np.zeros((Nt, Nx))
  U[0, :] = u0.copy()

  # Step 1: First time step uses Upwind to initialize U[1, :]
#   for j in range(1, Nx - 1):
#     U[1, j] = U[0, j] - alpha * (U[0, j] - U[0, j - 1])
#   U[1, 0] = 0.0
#   U[1, -1] = U[1, -2]

  # Remaining time steps use 3-level Leapfrog
  for i in range(1, Nt - 1):
    for j in range(1, Nx - 1):
      U[i + 1, j] = U[i - 1, j] - alpha * (U[i, j + 1] - U[i, j - 1])
    U[i + 1, 0] = 0.0
    U[i + 1, -1] = U[i + 1, -2]
  return U



# ==========================================
# 3. SOLVE ALL METHODS
# ==========================================
solutions = {
    "Upwind": solve_upwind(u0, alpha, Nt, Nx),
    "Lax-Friedrichs": solve_lax_friedrichs(u0, alpha, Nt, Nx),
    "Lax-Wendroff": solve_lax_wendroff(u0, alpha, Nt, Nx),
    "Leapfrog": solve_leapfrog(u0, alpha, Nt, Nx),
   
    # "FTCS": solve_ftcs(u0, alpha, Nt, Nx)  # Unstable, omitted from plot for clarity
}

# ==========================================
# 4. ANIMATED COMPARISON PLOT
# ==========================================
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_xlim(0, x[-1])
ax.set_ylim(-0.2, 1.3)
ax.set_xlabel("x (Position)")
ax.set_ylabel("U (Amplitude)")
ax.grid(True)

# Create line plot objects for each method
lines = {}
colors = {
    "Upwind": "blue",
    "Lax-Friedrichs": "green",
    "Lax-Wendroff": "red",
    "Leapfrog": "purple",
   
}

for method_name in solutions.keys():
  (line,) = ax.plot([], [], label=method_name, color=colors[method_name], lw=2)
  lines[method_name] = line

ax.legend(loc="upper right")


def init():
  for line in lines.values():
    line.set_data([], [])
  return list(lines.values())


def update(frame):
  for method_name, U_mat in solutions.items():
    lines[method_name].set_data(x, U_mat[frame, :])
  ax.set_title(
      f"Comparison of Advection Schemes (CFL = {alpha:.2f}) - Step:"
      f" {frame}/{Nt-1}"
  )
  return list(lines.values())


ani = FuncAnimation(
    fig, update, frames=Nt, init_func=init, interval=50, blit=True
)
# --------------------------------------------------
# INSERT BUTTON 
# --------------------------------------------------
ax_button = plt.axes([0.45, 0.02, 0.12, 0.05])
btn_pause = Button(ax_button, 'Pause / Play')

# Variable to track play/pause state
is_paused = False


# Function triggered when button is clicked
def toggle_pause(event):
  global is_paused
  if is_paused:
    ani.resume()
    is_paused = False
  else:
    ani.pause()
    is_paused = True


# Connect the button click event to the toggle function
btn_pause.on_clicked(toggle_pause)

plt.show()