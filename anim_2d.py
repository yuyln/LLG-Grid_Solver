import matplotlib.pyplot as plt;
import matplotlib.animation as ani;
import numpy as np

def cria(size):
    fig = plt.figure(figsize=size);
    ax = fig.add_subplot();
    return fig, ax;

def le_arq(path):
    with open(path, "r") as arq:
        linhas = arq.readlines();
        colunas = [[] for i in linhas[0].split(" \t ")[:-1]];
        for linha in linhas:
            for i, ele in enumerate(linha.split(" \t ")[:-1]):
                colunas[i].append(float(ele));
        for item in colunas:
            item = np.array(item)
        
    return np.array(colunas);

def pegavec(i):
        global colunas;
        global pos;
        x, y = zip(*[pos[j: j+2, i] for j in range(0, len(colunas), 3)]);
        u, v = zip(*[colunas[j: j+2, i] for j in range(0, len(colunas), 3)]);
        return x, y, u, v;
def pegacor(i):
    global colunas;
    cor = [colunas[j+2, i] for j in range(0, len(colunas), 3)];
    return np.array(cor)

def pegaheat(i):
    global max_
    return pegacor(i).reshape(max_, -1)

def anima(i):
        global b;
        global b_;
        b.remove();
        b_.remove()
        b = ax.quiver(*pegavec(i), linewidth=5)
        b_ = ax.imshow(pegaheat(i), cmap='bwr')

def faz(tempo):
        fps = int(len(colunas[0]) / tempo);
        a = ani.FuncAnimation(fig, anima, frames=len(colunas[0]));
        a.save("out.mp4", fps=fps, dpi=200);

colunas = le_arq("out.dat")
pos = le_arq("pos.dat")

fig, ax = cria([10, 10])
max_ = int((len(colunas) / 3) ** (1/2))
b = ax.quiver(*pegavec(100), linewidth=5)
b_ = ax.imshow(pegaheat(100), cmap='bwr')
b_.figure.colorbar(b_, ax=ax)
ax.set_xlim([-1, max_])
ax.set_ylim([-1, max_])
#ax.set_zlim([-1, 1])
#plt.show()
faz(20)
