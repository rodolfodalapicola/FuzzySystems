import matplotlib.pyplot as plt
import numpy as np

DISC_POINTS = 1000
FUZZY_LIMIT = 14


def idx(x):
    return int(x * DISC_POINTS / FUZZY_LIMIT)


def nmb(x):
    return x * FUZZY_LIMIT / DISC_POINTS


def plotAll(space, conjs):
    fig =  plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    if (len(conjs) == DISC_POINTS):
        plt.plot(space, conjs)
    else:
        for conj in conjs:
            plt.plot(space, conj, '*')
    ax.set_ylim([0, 1.1])
    plt.show()


def getSpaceVec():
    return np.array([i * FUZZY_LIMIT / DISC_POINTS for i in range(DISC_POINTS)])


def makeTriangular(linear_start, linear_top):
    a = (linear_top - linear_start) / 1
    x = np.zeros(DISC_POINTS)
    for i in range(idx(linear_top - linear_start)):
        if(0 <= i+idx(linear_start) < DISC_POINTS):
            x[i+idx(linear_start)] = nmb(i) / a
        if(0 <= i+idx(linear_top) < DISC_POINTS):
            x[i+idx(linear_top)] = 1 - nmb(i) / a
    return x


def makeTrap(linear_start, const_start, const_finish):
    a = (const_start - linear_start) / 1
    x = np.zeros(DISC_POINTS)
    for i in range(idx(const_start - linear_start)):
        if nmb(i) + linear_start >= 0:
            x[i+idx(linear_start)] = nmb(i) / a
        if nmb(i) + const_finish < FUZZY_LIMIT:
            x[i+idx(const_finish)] = 1 - nmb(i) / a
    for i in range(idx(const_finish - const_start)):
        if nmb(i) + const_start < FUZZY_LIMIT:
            x[i + idx(const_start)] = 1
    return x


def active(y, x):
    if y[idx(x)] != 0:
        return True
    return False


def getActives(x, conjs):
    active_fuzzy = []
    for conj in conjs:
        if active(conj, x):
            active_fuzzy.append(conj)
    return active_fuzzy


def getMembership(x, y):
    return y[idx(x)]


def crisp(conj, membership):
    crisped = conj
    for i in range(len(crisped)):
        if (crisped[i] < membership):
            crisped[i] = 0
    return crisped


def union(conjs):
    result = np.zeros(len(conjs[0]))
    for i in range(len(result)):
        arr = [0 for i in range(len(conjs))]
        for j in range(len(conjs)):
            arr[j] = conjs[j][i]
        result[i] = np.amax(arr)
    return result


def sumUnion(conjs):
    result = np.zeros(len(conjs[0]))
    for i in range(len(result)):
        arr = conjs[0][i]
        for j in range(1, len(conjs)):
            arr = arr + conjs[j][i] - arr * conjs[j][i]
        result[i] = arr
        if result[i] > 1:
            result[i] = 1
    return result


def activeIntersec(conjs):
    result = np.zeros(len(conjs[0]))
    for i in range(len(result)):
        act = getActives(nmb(i), conjs)
        if (len(act) == 0):
            result[i] = 0
            continue
        arr = [0 for i in range(len(act))]
        for j in range(len(act)):
            arr[j] = act[j][i]
        result[i] = np.amin(arr)
    return result


def intersec(conjs):
    result = np.zeros(len(conjs[0]))
    for i in range(len(result)):
        arr = [0 for i in range(len(conjs))]
        for j in range(len(conjs)):
            arr[j] = conjs[j][i]
        result[i] = np.amin(arr)
    return result


def prodIntersec(conjs):
    result = np.zeros(len(conjs[0]))
    for i in range(len(result)):
        arr = conjs[0][i]
        for j in range(1, len(conjs)):
            arr = arr * conjs[j][i]
        result[i] = arr
    return result


def complement(conj):
    result = np.zeros(len(conj))
    for i in range(len(result)):
        result[i] = 1 - conj[i]
    return result


def max_min(rxy, ryz):
    ros = np.zeros([len(rxy), len(ryz[0])])
    for x in range(len(rxy)):
        for z in range(len(ryz[0])):
            y_vec = np.zeros(len(ryz))
            for y in range(len(y_vec)):
                y_vec[y] = min(rxy[x][y], ryz[y][z])
            ros[x][z] = np.amax(y_vec)
    return ros


# ----------------------------------------------------------------------------------------------------------------------
# DAQUI PARA BAIXO ESTAO OS REQUERIMENTOS DO EPC 2
# ----------------------------------------------------------------------------------------------------------------------
def singleton(conj, x):
    output = np.zeros(len(conj))
    output[idx(x)] = 1
    return output


def mamdani(conj_a, conj_b):
    output = np.minimum(np.expand_dims(conj_a, axis=1), np.expand_dims(conj_b, axis=0))
    return output

def minimum(conj_a, conj_b):
    output = np.minimum(conj_a, conj_b)
    return output

def mamdani2(conj_a, conj_b):
    output = np.minimum(conj_b, conj_a)
    return output


def zadeh(conj_a, conj_b):
    output = np.minimum(np.expand_dims(conj_a, axis=1), np.expand_dims(conj_b, axis=0))
    max_comp = 1 - conj_a
    output = np.maximum(output, np.expand_dims(max_comp, axis=1))
    return output


def larsen(conj_a, conj_b):
    output = np.outer(conj_a, conj_b)
    return output

# ----------------------------------------------------------------------------------------------------------------------
# DAQUI PARA BAIXO ESTAO OS REQUERIMENTOS DO EPC 7 e do projeto
# ----------------------------------------------------------------------------------------------------------------------
def getAlfas(fuzzy_a, fuzzy_b, alfa):
    a1 = 0
    b1 = 0
    a2 = 0
    b2 = 0
    for x in range(len(fuzzy_a)):
        if (fuzzy_a[x] > alfa):
            a1 = x
            break
    for x in range(a1+1, len(fuzzy_a)):
        if (fuzzy_a[x] < alfa):
            b1 = x
            break
    for x in range(len(fuzzy_b)):
        if (fuzzy_b[x] > alfa):
            a2 = x
            break
    for x in range(a2+1, len(fuzzy_b)):
        if (fuzzy_b[x] < alfa):
            b2 = x
            break
    return a1, b1, a2, b2

def sumFuzzy(fuzzy_a, fuzzy_b, alfa_num):
    step = 1/alfa_num
    output = np.zeros(len(fuzzy_a))
    for i in range(alfa_num):
        a1, b1, a2, b2 = getAlfas(fuzzy_a, fuzzy_b, i*step)
        if a1+a2 < len(fuzzy_a):
            output[a1+a2] = i * step
        if b1+b2 < len(fuzzy_a):
            output[b1+b2] = i * step
    return output

def subFuzzy(fuzzy_a, fuzzy_b, alfa_num):
    step = 1/alfa_num
    output = np.zeros(len(fuzzy_a))
    over_0 = 0
    total = 0
    for i in range(alfa_num):
        alfa_corte = i * step
        a1, b1, a2, b2 = getAlfas(fuzzy_a, fuzzy_b, alfa_corte)
        if a1-b2 >= 0:
            output[a1-b2] = alfa_corte
        if b1-a2 >= 0:
            output[b1-a2] = alfa_corte

    med_area = 0
    for i in range(len(output)):
        med_area += output[i]

    med_area = med_area / 2
    soma = 0
    indice_goal = 0
    for i in range(len(output)):
        soma += output[i]
        if (soma > med_area):
            indice_goal = i
            break

    return output, indice_goal

def subConst(fuzzy, cnst):
    output = np.zeros(len(fuzzy))

    med_area = 0
    for i in range(cnst, len(output)):
        med_area += fuzzy[i]

    med_area = med_area / 2

    soma = 0
    indice_goal = 0
    for i in range(cnst, len(output)):
        soma += fuzzy[i]
        if (soma > med_area):
            indice_goal = i
            break
    output = output + getMembership(indice_goal, fuzzy)
    return output, indice_goal

def minFuzzy(fuzzy_a, fuzzy_b, alfa_num):
    step = 1/alfa_num
    output = np.zeros(len(fuzzy_a))
    over_0 = 0
    total = 0
    for i in range(alfa_num):
        alfa_corte = i * step
        a1, b1, a2, b2 = getAlfas(fuzzy_a, fuzzy_b, alfa_corte)
        output[min(a1, a2)] = alfa_corte
        output[min(b1, b2)] = alfa_corte

    return output

def maxFuzzy(fuzzy_a, fuzzy_b, alfa_num):
    step = 1/alfa_num
    output = np.zeros(len(fuzzy_a))
    over_0 = 0
    total = 0
    for i in range(alfa_num):
        alfa_corte = i * step
        a1, b1, a2, b2 = getAlfas(fuzzy_a, fuzzy_b, alfa_corte)
        output[max(a1, a2)] = alfa_corte
        output[max(b1, b2)] = alfa_corte

    return output

def prodFuzzy(fuzzy_a, fuzzy_b, alfa_num):
    step = 1/alfa_num
    output = np.zeros(len(fuzzy_a))
    for i in range(alfa_num):
        a1, b1, a2, b2 = getAlfas(fuzzy_a, fuzzy_b, i*step)
        if idx(nmb(a1)*nmb(a2)) < len(fuzzy_a):
            output[idx(nmb(a1)*nmb(a2))] = i * step
        if idx(nmb(b1)*nmb(b2)) < len(fuzzy_a):
            output[idx(nmb(b1)*nmb(b2))] = i * step
    return output

def divFuzzy(fuzzy_a, fuzzy_b, alfa_num):
    step = 1/alfa_num
    output = np.zeros(len(fuzzy_a))
    for i in range(alfa_num):
        a1, b1, a2, b2 = getAlfas(fuzzy_a, fuzzy_b, i*step)
        b2 += 1
        if idx(nmb(a1)/nmb(b2)) < len(fuzzy_a):
            output[idx(a1/b2)] = i * step
        if idx(nmb(b2)/nmb(a1)) < len(fuzzy_a):
            output[idx(nmb(b2)/nmb(a1))] = i * step
    return output

def defuzzify(fuzzy):
    return max(fuzzy) 
