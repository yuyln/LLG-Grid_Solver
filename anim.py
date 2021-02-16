import matplotlib.pyplot as plt;
import matplotlib.animation as ani;
import numpy as np

def cria(size):
    fig = plt.figure(figsize=size);
    ax = fig.add_subplot(111, projection='3d');
    return fig, ax;

def le_arq(path):
    with open(path, "r") as arq:
        linhas = arq.readlines();
        colunas = [[] for i in linhas[0].split(" \t ")[:-1]];
        for linha in linhas:
            for i, ele in enumerate(linha.split(" \t ")[:-1]):
                colunas[i].append(float(ele));
    return np.array(colunas);

def pegavec(i):
        global colunas;
        global pos;
        x, y, z = zip(*[pos[j: j+3, i] for j in range(0, len(colunas), 3)]);
        u, v, w = zip(*[colunas[j: j+3, i] for j in range(0, len(colunas), 3)]);
        return x, y, z, u, v, w;

def anima(i):
        global b;
        b.remove();
        b = ax.quiver(*pegavec(i), linewidth=5);

def faz(tempo):
        fps = int(len(colunas[0]) / tempo);
        a = ani.FuncAnimation(fig, anima, frames=len(colunas[0]), interval=1);
        a.save("out.mp4", fps=fps, dpi=200);

colunas = le_arq("out.dat")
pos = le_arq("pos.dat")

fig, ax = cria([10, 10])
b = ax.quiver(*pegavec(0), linewidth=5)
ax.set_xlim([-1, 11])
ax.set_ylim([-1, 11])
ax.set_zlim([-1, 11])
#plt.show()
faz(10)
