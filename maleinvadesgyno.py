import numpy as np
import math
def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper

import matplotlib.pyplot as plt
import pandas as pd

cs = 0.0 #cost of the cytoplasmic sterility
bigscreen = []
timetoequilibrium = []
for pollen_threshold in np.around([1], decimals=4):
    for csm in np.around([0.8], decimals=4):
        for alpha in np.around(np.linspace(0.5, 16, 55), decimals=4):  # np.linspace(0.5,60,200)
            for s in np.around([0.4], decimals=4):  # np.linspace(0,0.95,40) 0.3
                for d in np.around([0.1], decimals=4):  # np.linspace(0,0.95,40) 0.1
                    for g in np.around(np.linspace(0.7, 10, 45), decimals=4):  # np.linspace(0.4,1.7,40)
                        Aa = 0  # frequency of male initial
                        AAc = 0.000001
                        AA = 1 - AAc  # frequency of hermaphrodites initial
                        Aac = 0

                        i = 0
                        stored_values = []
                        stored_values.append([i, Aa, AA, Aac, AAc, s, d])
                        Condition = False
                        while Condition == False:
                            i = i + 1

                            L = min(1, AA / 1)
                            L = 1  # without pollen limitation
                            pollen_a = 0
                            pollen_A = 1
                            AA_number = (L * AA * (1 - s) * pollen_A
                                         + AA * s * (1 - d) ) # fertilized ovules not suffer from polim
                                         #+ AA * (1 - L) * (1 - s) * (1 - d)) #re fertilize the pollim ovules
                            Aa_number = L * AA * (1 - s) * pollen_a
                            AAc_number = L * AAc * g * pollen_A
                            Aac_number = L * AAc * g * pollen_a
                            total = AA_number + Aa_number + AAc_number + Aac_number
                            if total == 0:
                                AA = 0
                                Aa = 0
                                Aac = 0
                                AAc = 0
                            else:
                                Aa = Aa_number / total
                                AA = AA_number / total
                                AAc = AAc_number / total
                                Aac = Aac_number / total

                            stored_values.append([i, Aa, AA, Aac, AAc, s, d])
                            if AAc == 0 or (
                                    stored_values[-2][2] - stored_values[-1][2] < 1E-10 and stored_values[-1][2] -
                                    stored_values[-2][
                                        2] < 1E-10):
                                Condition = True

                        AAc = min((1 - AA), 0.99)
                        Aac = 0
                        Aa = 0.000001
                        AA = 1 - Aa - AAc

                        j = 0
                        stored_values = []
                        stored_values.append([j, Aa, AA, Aac, AAc, s, d])
                        Condition = False
                        while Condition == False:
                            j = j + 1
                            total_pollen = (AA + Aa * alpha + Aac * alpha * (1 - cs) * (1 - csm))
                            if round(total_pollen, 4) == 0:
                                AA = 0
                                Aa = 0
                                Aac = 0
                                AAc = 0
                            else:
                                L = min(1, total_pollen / pollen_threshold)
                                L = 1 #without pollen limitation
                                pollen_a = (Aa * alpha * 0.5
                                            + Aac * alpha * (1 - cs) * (1 - csm) * 0.5) / total_pollen
                                pollen_A = 1 - pollen_a
                                AA_number = (L * AA * (1 - s) * pollen_A
                                             + AA * s * (1 - d) ) #fertilized ovules not suffer from polim
                                             #+ AA * (1 - L) * (1 - s) * (1 - d)) #re fertilize the pollim ovules
                                Aa_number = L * AA * (1 - s) * pollen_a
                                AAc_number = L * AAc * g * pollen_A
                                Aac_number = L * AAc * g * pollen_a
                                total = AA_number + Aa_number + AAc_number + Aac_number
                                if total == 0:
                                    AA = 0
                                    Aa = 0
                                    Aac = 0
                                    AAc = 0
                                else:
                                    Aa = Aa_number / total
                                    AA = AA_number / total
                                    AAc = AAc_number / total
                                    Aac = Aac_number / total
                            stored_values.append([j, Aa, AA, Aac, AAc, s, d])
                            if AAc == 0 or (j >1000  and stored_values[-20][2] - stored_values[-1][2] < 1E-10 and stored_values[-1][2] - stored_values[-20][
                                        2] < 1E-10):
                                Condition = True
                        bigscreen.append([truncate(Aa, 8),
                                          truncate(AA, 8),
                                          truncate(Aac, 8),
                                          truncate(AAc, 8),
                                          s, d, alpha, g, csm])


#### Plot the data
stored_values = pd.DataFrame(bigscreen)
stored_values.columns = ["Aa", "AA", "Aac", "AAc", "s", "d", "alpha", "g", "csm"]
stored_values.iloc[:, 0:4] = np.where(stored_values.iloc[:, 0:4] > 0.0, 1.0, 0.0)
stored_values["ID"] = stored_values.iloc[:, 0].astype(str) + stored_values.iloc[:, 1].astype(str) + stored_values.iloc[:, 2].astype(str) + stored_values.iloc[:, 3].astype(str)
stored_values.columns = ["Aa", "AA", "Aac", "AAc", "s", "d", "alpha", "g", "csm", "ID"]
systems = {"1.01.01.01.0": ["Trioecy", "purple"],
           "0.00.01.01.0": ["Dioecy", "blue"],
           "1.01.00.00.0": ["Androdioecy", "green"],
           "0.00.00.00.0": ["Extinction", "black"],
           "0.01.00.00.0": ["Hermaphroditism", "yellow"],
           "0.01.00.01.0": ["Gynodioecy", "orange"],
           "1.01.00.01.0": ["Trioecy", "purple"],
           "0.01.01.01.0": ["Trioecy", "purple"]}
systemscolor = {"Trioecy": "purple",
                "Dioecy": "blue",
                "Androdioecy": "green",
                "Extinction": "black",
                "Hermaphroditism": "yellow",
                "Gynodioecy": "orange",
                "cyto Trioecy": "pink",
                "perfect Trioecy": "cyan"}
systemscytinvade = {"Trioecy": "yes",
                    "Dioecy": "yes",
                    "Androdioecy": "no",
                    "Extinction": "yes",
                    "Hermaphroditism": "no",
                    "Gynodioecy": "yes",
                    "cyto Trioecy": "yes",
                    "perfect Trioecy": "yes"}
fateinvade = {"Trioecy": "polymorphic",
              "Gynodioecy": "polymorphic",
              "Dioecy": "fixed",
              "Androdioecy": "no",
              "Extinction": "no",
              "Hermaphroditism": "no",
              "cyto Trioecy": "polymorphic",
              "perfect Trioecy": "polymorphic"}


stored_values["system"] = [systems[i][0] for i in stored_values["ID"]]
stored_values["color"] = [systems[i][1] for i in stored_values["ID"]]
stored_values["invaded"] = [systemscytinvade[i] for i in stored_values["system"]]
stored_values["CMSfate"] = [fateinvade[i] for i in stored_values["system"]]


sub = stored_values
plt.figure(figsize = [6.13,5.03])
for sy in np.unique(stored_values["system"]):
    ix = np.where((sub["system"] == sy))
    plt.scatter(sub.loc[ix, "alpha"],sub.loc[ix, "g"],c=systemscolor[sy],alpha=0.5,label=sy)

plt.plot([2 * (1 - s * d) / (1 - s)] * 100, np.linspace(0.7, 10, 100), lw=5, c="black")
plt.plot([2 * (1 - s * d) / (1 - s)] * 100, np.linspace(0.7, 3, 100), lw=5, c="black")
x = np.linspace(0.5,2*(1-s*d)/(1-s),100)
plt.plot(x, [1-s*d]*100, lw = 5, c = "red")
x = np.linspace(2*(1-s*d)/(1-s),16,100)
y = np.linspace(0.7,3,100)
Aaeq = (x*(1-s)-2*(1-s*d))/(2*(x-1)*(1-s*d))
polA = (1 - Aaeq + Aaeq * x * 0.5) / (1- Aaeq + Aaeq * x)
plt.plot(x, ((1-s)*polA+s*(1-d))/polA, lw = 5, c = "red")

x = np.linspace(0.5,16,100)
poly = 1-s+s*(1-d)/(0.25*x*(1-csm))
poly2 = [1+s-2*s*d]*100
poly3 = np.maximum(poly,poly2)
poly4 = (1-s)*(2*d-1)+(1-d)/(0.25*x*(1-csm))
poly5 = np.maximum(poly2,poly4)
plt.plot(x, poly2, lw = 5, c = "green")
plt.plot(x, poly3, lw = 5, c = "green")
plt.plot(x, poly5, lw = 5, c = "green")

Ycal = (x - 2)/(x-2*d)
Xcal = (Ycal*(1-d)*s)/(Ycal*(1+s)-1)
Yfin = Xcal + 1 - s
plt.plot(x, Yfin, lw = 3, c = "black")

plt.xlabel(r'Male fertility')
plt.ylabel("Female fertility")
plt.title(r's = {}, d = {}, ec = {}, with pollen limitation'.format(s, d, csm))
plt.ylim(0.585, 3.115)
plt.ylim(0.2, 10.45)