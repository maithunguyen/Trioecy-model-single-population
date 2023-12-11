import numpy as np
import math
def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper
import matplotlib.pyplot as plt
import pandas as pd

#s = 0.0 #proportion of self fertilized in hermaphrodites
#d = 0.1 #inbreeding depression
#csm = 0.8 #cost of the cytoplasmic sterility that makes male sterile (separated from cs)
#alpha = 12.52040
#g = 1.04290
cs = 0.0 #cost of the cytoplasmic sterility
bigscreen = []
timetoequilibrium = []
for csm in np.around([1], decimals=4): #cost of the cytoplasmic sterility on males
    for alpha in np.around(np.linspace(0.5,16,55), decimals=4): #male fertility advantage
        for s in np.around([0.4], decimals=4): #proportion of self fertilized in hermaphrodites
            for d in np.around([0.1], decimals=4): #inbreeding depression
                for g in np.around(np.linspace(0.7,1.4,45), decimals=4): #female fertility advantage
                    Aa = 0.0001  # frequency of male initial
                    AA = 1 - Aa  # frequency of hermaphrodites initial
                    AAc = 0
                    Aac = 0

                    i = 0
                    stored_values = []
                    stored_values.append([i, Aa, AA, Aac, AAc, s, d])
                    Condition = False
                    while Condition == False:
                        i = i + 1
                        pollen_a = (Aa * alpha * 0.5) / (AA + Aa * alpha)
                        pollen_A = 1 - pollen_a
                        AA_new = AA * s * (1 - d) + AA * (1 - s) * pollen_A
                        Aa_new = AA * (1 - s) * pollen_a
                        AA = AA_new / (AA_new + Aa_new)
                        Aa = 1 - AA
                        stored_values.append([i, Aa, AA, Aac, AAc, s, d])
                        if Aa == 0 or (stored_values[-2][2] - stored_values[-1][2] < 1E-10 and stored_values[-1][2] - stored_values[-2][
                        2] < 1E-10):
                            Condition = True

                    Aa = truncate(Aa, 8)
                    Aac = 0
                    AAc = 0.000001
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
                            pollen_a = (Aa * alpha * 0.5
                                    + Aac * alpha * (1 - cs) * (1 - csm) * 0.5) / total_pollen
                            pollen_A = 1 - pollen_a
                            AA_number = (AA * (1 - s) * pollen_A
                                     + AA * s * (1 - d))
                            Aa_number = AA * (1 - s) * pollen_a
                            AAc_number = AAc * g * pollen_A
                            Aac_number = AAc * g * pollen_a
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
                        if AAc == 0 or (j > 20 and stored_values[-2][2] - stored_values[-1][2] < 1E-10 and stored_values[-1][2] - stored_values[-2][
                                2] < 1E-10):
                                Condition = True
                    timetoequilibrium.append([i, j, s, d, alpha, g, csm])
                    bigscreen.append([truncate(Aa, 8),
                                    truncate(AA, 8),
                                    truncate(Aac, 8),
                                    truncate(AAc, 8),
                                    s, d, alpha, g, csm])

# Plot the data
stored_values = pd.DataFrame(bigscreen)
stored_values.columns = ["Aa", "AA", "Aac", "AAc", "s", "d", "alpha", "g", "csm"]
# stored_values. to_pickle("bigscreen3.pkl")
# stored_values = pd. read_pickle("bigscreen3.pkl")
stored_values.iloc[:, 0:4] = np.where(stored_values.iloc[:, 0:4] > 0.000001, 1.0, 0.0)
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

#output = pd.DataFrame(bigscreen)
#output.columns = ["Aa", "AA", "Aac", "AAc", "s", "d", "alpha", "g", "csm"]
#output["system"] = stored_values["system"]
#output["invaded"] = stored_values["invaded"]
#output["CMSfate"] = stored_values["CMSfate"]

#output.to_csv(r'alpha_g_s03d01ec0.csv', index = False, header=True)
sub = stored_values
plt.figure(figsize = [6.13,5.03])
for sy in np.unique(stored_values["system"]):
    ix = np.where((sub["system"] == sy))
    plt.scatter(sub.loc[ix, "alpha"],sub.loc[ix, "g"],c=systemscolor[sy],alpha=0.5,label=sy)

plt.plot([2 * (1 - s * d) / (1 - s)] * 100, np.linspace(0.7, 1.4, 100), lw=5, c="black")
x = np.linspace(0.5,2*(1-s*d)/(1-s),100)
plt.plot(x, [1-s*d]*100, lw = 5, c = "red")
x = np.linspace(2*(1-s*d)/(1-s),16,100)
Aaeq = (x*(1-s)-2*(1-s*d))/(2*(x-1)*(1-s*d))
polA = (1 - Aaeq + Aaeq * x * 0.5) / (1- Aaeq + Aaeq * x)
plt.plot(x, ((1-s)*polA+s*(1-d))/polA, lw = 5, c = "red")
plt.plot(x, [1+s-2*s*d]*100, lw = 5, c = "green")


plt.annotate('CMS fixed', xy=(8, 1.18), xycoords='data',bbox=dict(boxstyle="round", fc="0.8"))
#plt.annotate('CMS\nfixed', xy=(0.5, 1.2), xycoords='data',bbox=dict(boxstyle="round", fc="0.8"))
#plt.annotate('CMS not invaded', xy=(8, 0.8), xycoords='data',bbox=dict(boxstyle="round", fc="0.8"))
#plt.annotate('CMS\nnot\ninvaded', xy=(0.5, 0.8), xycoords='data',bbox=dict(boxstyle="round", fc="0.8"))
#plt.annotate('CMS invaded', xy=(8, 1.18), xycoords='data',bbox=dict(boxstyle="round", fc="0.8"))
plt.annotate('CMS\nfixed', xy=(0.5, 1), xycoords='data',bbox=dict(boxstyle="round", fc="0.8"))
plt.annotate('CMS fixed', xy=(8, 1.3), xycoords='data',bbox=dict(boxstyle="round", fc="0.8"))
#plt.annotate('Worse', xy=(8, 1.12), xytext=(8, 1.18)),arrowprops=dict(arrowstyle='<-'), bbox=dict(boxstyle="round", fc="0.8"))
plt.annotate('CMS\nnot\ninvaded', xy=(0.5, 0.8), xycoords='data',bbox=dict(boxstyle="round", fc="0.8"))
plt.annotate('CMS not invaded', xy=(8, 1), xycoords='data',bbox=dict(boxstyle="round", fc="0.8"))
plt.legend(loc='center right')
plt.xlabel(r'$\alpha$')
plt.ylabel("g")
plt.title(r'Outcome of CMS mutation invasion, s = {}, d = {}, ec = {}'.format(s, d, csm))

#### Plot lower female
stored_values = pd.DataFrame(bigscreen)
stored_values.columns = ["Aa", "AA", "Aac", "AAc", "s", "d", "alpha", "g", "csm"]
stored_values.loc[(stored_values.AAc < 0.000001), "AAc"] = 0.0
stored_values.loc[(round(stored_values.AAc,8) >0.999999), "AAc"] = 1.0
stored_values.loc[(stored_values.AAc > 0.0) & (stored_values.AAc <= 0.2), "AAc"] = 2
stored_values.loc[(stored_values.AAc > 0.2) & (stored_values.AAc <= 0.4), "AAc"] = 4
stored_values.loc[(stored_values.AAc > 0.4) & (round(stored_values.AAc,8) <= 0.5), "AAc"] = 5
stored_values.loc[(round(stored_values.AAc,8) > 0.5) & (stored_values.AAc <= 0.6), "AAc"] = 6
stored_values.loc[(stored_values.AAc > 0.6) & (stored_values.AAc <= 0.8), "AAc"] = 8
stored_values.loc[(stored_values.AAc > 0.8) & (round(stored_values.AAc,8) < 0.999999), "AAc"] = 9


stored_values["ID"] = stored_values.iloc[:, 3].astype(str)
stored_values.columns = ["Aa", "AA", "Aac", "AAc", "s", "d", "alpha", "g", "csm", "ID"]
systems = {"0.0": ["0", "grey"],
            "2.0": ["0 - 0.2", "pink"],
            "4.0": ["0.2 - 0.4", "red"],
            "5.0": ["0.4 - 0.5", "yellow"],
            "6.0": ["0.5 - 0.6", "green"],
            "8.0": ["0.6 - 0.8", "orange"],
           "9.0": ["0.8 - <1", "blue"],
           "1.0": ["1", "black"]}
systemscolor = {"0": "grey",
                "0 - 0.2": "pink",
                "0.2 - 0.4": "red",
                "0.4 - 0.5": "yellow",
                "0.5 - 0.6": "green",
                "0.6 - 0.8": "orange",
                "0.8 - <1": "blue",
                "1": "black"}


stored_values["system"] = [systems[i][0] for i in stored_values["ID"]]
stored_values["color"] = [systems[i][1] for i in stored_values["ID"]]

sub = stored_values
for sy in np.unique(stored_values["system"]):
    ix = np.where((sub["system"] == sy))
    plt.scatter(sub.loc[ix, "alpha"],sub.loc[ix, "g"],c=systemscolor[sy],alpha=0.5,label=sy)
x = np.linspace(0.5,2*(1-s*d)/(1-s),100)
plt.plot(x, [1-s*d]*100, lw = 3, c = "red")
x = np.linspace(2*(1-s*d)/(1-s),20,100)
Aaeq = (x*(1-s)-2*(1-s*d))/(2*(x-1)*(1-s*d))
polA = (1 - Aaeq + Aaeq * x * 0.5) / (1- Aaeq + Aaeq * x)
plt.plot(x, ((1-s)*polA+s*(1-d))/polA, lw = 3, c = "red")
poly = 1-s+s*(1-d)/(0.25*x*(1-csm))
poly2 = [1+s-2*s*d]*100
poly3 = np.maximum(poly,poly2)
plt.plot(x, poly2, lw = 3, c = "green")
#plt.plot(x, poly2, lw = 3, c = "green")
plt.plot([2*(1-s*d)/(1-s)]*100, np.linspace(0.7,1.5,100), lw = 3, c = "pink")
plt.legend(loc='lower right')
plt.xlabel(r'Male fertility')
plt.ylabel("Female fertility")
plt.title(r'Female frequency, s = {}, d = {}, ec = {}'.format(s, d, csm))

#### Plot CMS
stored_values = pd.DataFrame(bigscreen)
stored_values.columns = ["Aa", "AA", "Aac", "AAc", "s", "d", "alpha", "g", "csm"]
stored_values["CMS"] = stored_values.AAc + stored_values.Aac
stored_values.loc[(stored_values.CMS == 0), "CMS"] = 0.0
stored_values.loc[(round(stored_values.CMS,8) >0.999999), "CMS"] = 1.0
stored_values.loc[(stored_values.CMS > 0.0) & (stored_values.CMS <= 0.2), "CMS"] = 2
stored_values.loc[(stored_values.CMS > 0.2) & (stored_values.CMS <= 0.4), "CMS"] = 4
stored_values.loc[(stored_values.CMS > 0.4) & (stored_values.CMS <= 0.6), "CMS"] = 6
stored_values.loc[(stored_values.CMS > 0.6) & (stored_values.CMS <= 0.8), "CMS"] = 8
stored_values.loc[(stored_values.CMS > 0.8) & (round(stored_values.CMS,8) < 0.999999), "CMS"] = 9
stored_values["ID"] = stored_values.iloc[:, 9].astype(str)
stored_values.columns = ["Aa", "AA", "Aac", "AAc", "s", "d", "alpha", "g", "csm", "CMS", "ID"]
systems = {"0.0": ["0", "grey"],
            "2.0": ["0 - 0.2", "pink"],
            "4.0": ["0.2 - 0.4", "red"],
            "6.0": ["0.4 - 0.6", "yellow"],
            "8.0": ["0.6 - 0.8", "orange"],
           "9.0": ["0.8 - <1", "blue"],
           "1.0": ["1", "black"]}
systemscolor = {"0": "grey",
                "0 - 0.2": "pink",
                "0.2 - 0.4": "red",
                "0.4 - 0.6": "yellow",
                "0.6 - 0.8": "orange",
                "0.8 - <1": "blue",
                "1": "black"}


stored_values["system"] = [systems[i][0] for i in stored_values["ID"]]
stored_values["color"] = [systems[i][1] for i in stored_values["ID"]]

sub = stored_values
for sy in np.unique(stored_values["system"]):
    ix = np.where((sub["system"] == sy))
    plt.scatter(sub.loc[ix, "alpha"],sub.loc[ix, "g"],c=systemscolor[sy],alpha=0.5,label=sy)
x = np.linspace(0.5,2*(1-s*d)/(1-s),100)
plt.plot(x, [1-s*d]*100, lw = 3, c = "red")
x = np.linspace(2*(1-s*d)/(1-s),20,100)
Aaeq = (x*(1-s)-2*(1-s*d))/(2*(x-1)*(1-s*d))
polA = (1 - Aaeq + Aaeq * x * 0.5) / (1- Aaeq + Aaeq * x)
plt.plot(x, ((1-s)*polA+s*(1-d))/polA, lw = 3, c = "red")
poly = 1-s+s*(1-d)/(0.25*x*(1-csm))
poly2 = [1+s-2*s*d]*100
poly3 = np.maximum(poly,poly2)
#plt.plot(x, poly3, lw = 3, c = "green")
plt.plot(x, poly2, lw = 3, c = "green")
plt.plot([2*(1-s*d)/(1-s)]*100, np.linspace(0.7,3,100), lw = 3, c = "pink")