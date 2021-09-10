import Plots as plt

pos_0 = [0, 0]
vel_0 = [10, 100]
N = 1000
dt = 0.02
g = 9.81

pos = zeros(2, N)
vel = zeros(2, N)
acc = zeros(2, N)

pos[:, 1] = pos_0
vel[:, 1] = vel_0;

grav(pos) = [0, -g]

for i in 1:N-1
    acc[:, i] = grav(pos[:, i])
    vel[:, i + 1] = vel[:, i] + acc[:, i] * dt
    pos[:, i + 1] = pos[:, i] + vel[:, i] * dt
end

plt.display(plt.plot(pos[2, :], line = 3))