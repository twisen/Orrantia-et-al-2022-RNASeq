import numpy as np;
from math import log, e;
from scipy.stats import trim_mean, ttest_ind;
from scipy.stats import pearsonr, spearmanr;

import matplotlib
matplotlib.use('Agg')


import matplotlib.pyplot as plt;
from statsmodels.sandbox.stats.multicomp import multipletests;

import os

directorio = "Figures/"

try:
	os.stat(directorio)
except:
	os.mkdir(directorio)

LIMITE_EXPR = 0.5;


NOMBRES = [];

CON_NOMBRE_NORMAL = True;
tipoNombre = "GENE_ID";
if CON_NOMBRE_NORMAL:
	tipoNombre = "GENE_ID"
	nombreGen = {};
	GTF = "../Clean_GTF/Homo_sapiens.GRCh38.91.gtf";
	f = open(GTF, "r");
	while True:
		linea = f.readline();
		if not linea: break;
		l = linea.split("\t");
		if len(l) < 7: continue;
		
		if tipoNombre == "GENE_ID":
			if l[2] != "gene": continue;
			datosGen = l[8].split("; ");
			ID = datosGen[0].split(" ")[1][1:-1];
			nombre = datosGen[2].split(" ")[1][1:-1];
		elif tipoNombre == "TRANSCRIPT_ID":
			if l[2] != "transcript": continue;
			datosGen = l[8].split("; ");
			ID = datosGen[2].split(" ")[1][1:-1];
			nombre = datosGen[4].split(" ")[1][1:-1];
		nombreGen[ID] = nombre;
	f.close();


fichero = "../Raw_Data/Raw_Data.txt";
"Reading file: " + fichero;
f = open(fichero, "r");
linea = "\t"+f.readline();
#~ linea = f.readline();
DATOS = {};
l = linea.split("\n")[0].split("\r")[0].split("\t");
transcritos = [];
nombres = []
for n in l[1:]:
	nombre = n.split("/")[-1].split(".")[0]
	if len(nombre) < 1: continue;
	nombres.append(nombre);
	DATOS[nombre] = {};
j = 0;
while True:
	linea = f.readline()
	if not linea: break;
	l = linea.split("\t");
	if len(l) < 2: continue;
	
	#~ transcrito = l[0][1:-1];
	transcrito = l[0];
	aux = transcrito.split("\"")
	if len(aux) > 1: transcrito = aux[1]
	transcritos.append(transcrito);
	for i in range(1, len(l)):
		DATOS[nombres[i-1]][transcrito] = int(float(l[i]));#float(l[i]);#int(float(l[i]));
	j+=1;
	if j%10000 == 0: print j;
f.close();
print j;
print "EOF";

print "NAMES: " + str(nombres);

TrancritosPorEXP = {};
EXPPorTranscrito = {};
for t in transcritos: EXPPorTranscrito[t] = 0;
iteraciones = 0;
for n in nombres:
	TrancritosPorEXP[n] = 0;
	for t in EXPPorTranscrito:
		if DATOS[n][t] > 0:
			EXPPorTranscrito[t] += 1;
			TrancritosPorEXP[n] += 1;
transcritosFiltrados = [];
LIM_EXP = 8;
for t in EXPPorTranscrito:
	if EXPPorTranscrito[t] >= (LIM_EXP): transcritosFiltrados.append(t);

print "Detected transcripts: " + str(len(transcritosFiltrados));

transcritos = transcritosFiltrados; ######Eliminamos los transcritos no detectados

TABLA = {};
TABLA["Nombres"] = transcritos;
for n in nombres:
	TABLA[n] = [];
	for t in transcritos: TABLA[n].append(DATOS[n][t]);

NOMBRES += nombres;

################Normalizaciones

FactoresNorm = {}

RLE = {};
RLE["Nombres"] = transcritos;
for n in nombres: RLE[n] = [];
for i in range(0, len(transcritos)):
	aux = [];
	for n in nombres: aux.append(TABLA[n][i]);
	aux = log(np.median(aux)+1, 2);
	for n in nombres: RLE[n].append(log(TABLA[n][i]+1, 2)-aux);

aux = [];
for n in nombres: aux.append(RLE[n]);
plt.title("RLE Raw");
plt.boxplot(aux, sym='');
plt.ylim(-5, 5);
plt.xticks(range(1, len(nombres)+1),nombres, rotation=27);
plt.savefig("Figures/RLE_Raw.png");


aux = [];
for n in nombres: aux.append(TABLA[n]);
plt.figure();
plt.title("Expression Raw");
plt.boxplot(aux, sym='');
plt.xticks(range(1, len(nombres)+1),nombres, rotation=27);
plt.savefig("Figures/Expression_Raw.png");

print "NORMALIZING"
#TPM
NORM = "TPM"
FactoresNorm[NORM] = {}
for n in NOMBRES:
	factor = sum(TABLA[n]) / 1000000.0;
	FactoresNorm[NORM][n] = factor;
	for i in range(0, len(TABLA[n])): TABLA[n][i] /= factor;

fout = open("IdentifiedTranscripts.txt", "w");
fout.write("Name\tIdentified\n");
GenesIdentificados = {};
#Transcritos identificados
for n in NOMBRES:
	cuenta = 0;
	GenesIdentificados[n] = [];
	for i in range(0, len(TABLA[n])):
		if(TABLA[n][i] >= 1):
			cuenta +=1;
			GenesIdentificados[n].append(TABLA["Nombres"][i]);
	GenesIdentificados[n] = set(GenesIdentificados[n]);
	print "Identified in "+ n + " : " + str(cuenta);
	fout.write(n + "\t" + str(cuenta) + "\n");
fout.close();


#TMM
ESCALA=1;
NORM = "TMM"
FactoresNorm[NORM] = {};
for n in NOMBRES:
	factor = trim_mean(TABLA[n], 0.1)/ESCALA;
	FactoresNorm[NORM][n] = factor;
	for i in range(0, len(TABLA[n])): TABLA[n][i] /= factor;


RLE = {};
RLE["Nombres"] = transcritos;
for n in nombres: RLE[n] = [];
for i in range(0, len(transcritos)):
	aux = [];
	for n in nombres: aux.append(TABLA[n][i]);
	aux = log(np.median(aux)+1, 2);
	for n in nombres: RLE[n].append(log(TABLA[n][i]+1, 2)-aux);

aux = [];
for n in nombres: aux.append(RLE[n]);
plt.figure();
plt.title("RLE TMM");
plt.boxplot(aux, sym='');
plt.ylim(-2, 2);
plt.xticks(range(1, len(nombres)+1),nombres, rotation=27);
plt.savefig("Figures/RLE_TMM.png");


aux = [];
for n in nombres: aux.append(TABLA[n]);
plt.figure();
plt.title("Expression TMM");
plt.boxplot(aux, sym='');
plt.ylim(-5, 15*ESCALA);
plt.xticks(range(1, len(nombres)+1),nombres, rotation=27);
plt.savefig("Figures/Expression_TMM.png");


fout = open("FactoresNorm.txt", "w");
fout.write("Sample")
for norm in FactoresNorm: fout.write("\t"+norm);
fout.write("\n");
for n in nombres:
	fout.write(n);
	for norm in FactoresNorm: fout.write("\t"+str(FactoresNorm[norm][n]));
	fout.write("\n");
fout.close();


MT_GENES=[
"ENSG00000210144", "ENSG00000210077", "ENSG00000198938", "ENSG00000210049", "ENSG00000198695",
"ENSG00000210107", "ENSG00000210164", "ENSG00000210140", "ENSG00000198886", "ENSG00000210117", 
"ENSG00000212907", "ENSG00000210191", "ENSG00000210156", "ENSG00000210154", "ENSG00000211459", 
"ENSG00000210195", "ENSG00000210174", "ENSG00000210082", "ENSG00000198899", "ENSG00000198804", 
"ENSG00000198727", "ENSG00000198840", "ENSG00000209082", "ENSG00000198786", "ENSG00000210112", 
"ENSG00000210151", "ENSG00000210135", "ENSG00000210127", "ENSG00000210100", "ENSG00000198712", 
"ENSG00000198763", "ENSG00000210176", "ENSG00000210194", "ENSG00000210196", "ENSG00000228253", 
"ENSG00000198888", "ENSG00000210184"
];
PorcentajeMT = {};

RIBO_GENES=[
"ENSG00000252422", "ENSG00000201920", "ENSG00000252745", "ENSG00000252072", "ENSG00000252539", "ENSG00000252260", "ENSG00000253077", "ENSG00000252987", "ENSG00000222455", 
"ENSG00000251904", "ENSG00000201968", "ENSG00000201109", "ENSG00000212258", "ENSG00000223138", "ENSG00000252864", "ENSG00000212497", "ENSG00000252002", "ENSG00000222747", 
"ENSG00000202056", "ENSG00000222248", "ENSG00000252653", "ENSG00000223092", "ENSG00000252760", "ENSG00000202174", "ENSG00000222268", "ENSG00000201394", "ENSG00000212336", 
"ENSG00000252927", "ENSG00000223274", "ENSG00000199276", "ENSG00000252942", "ENSG00000223318", "ENSG00000223131", "ENSG00000200709", "ENSG00000199609", "ENSG00000199535", 
"ENSG00000202430", "ENSG00000201736", "ENSG00000252845", "ENSG00000272435", "ENSG00000199556", "ENSG00000253058", "ENSG00000200527", "ENSG00000253039", "ENSG00000238379", 
"ENSG00000252591", "ENSG00000252806", "ENSG00000252667", "ENSG00000276861", "ENSG00000201618", "ENSG00000199771", "ENSG00000212289", "ENSG00000252368", "ENSG00000252182", 
"ENSG00000252179", "ENSG00000200168", "ENSG00000252964", "ENSG00000252307", "ENSG00000200301", "ENSG00000251770", "ENSG00000252510", "ENSG00000199404", "ENSG00000200225", 
"ENSG00000252149", "ENSG00000252107", "ENSG00000222209", "ENSG00000276599", "ENSG00000201846", "ENSG00000200326", "ENSG00000222129", "ENSG00000201671", "ENSG00000223169", 
"ENSG00000212559", "ENSG00000212242", "ENSG00000252302", "ENSG00000199837", "ENSG00000239152", "ENSG00000252816", "ENSG00000252231", "ENSG00000222230", "ENSG00000222546", 
"ENSG00000251884", "ENSG00000223046", "ENSG00000207129", "ENSG00000202044", "ENSG00000252563", "ENSG00000212290", "ENSG00000201766", "ENSG00000200626", "ENSG00000222150", 
"ENSG00000201086", "ENSG00000251965", "ENSG00000222854", "ENSG00000238965", "ENSG00000274467", "ENSG00000200755", "ENSG00000222468", "ENSG00000200650", "ENSG00000222308", 
"ENSG00000222806", "ENSG00000252370", "ENSG00000212308", "ENSG00000252167", "ENSG00000251915", "ENSG00000252428", "ENSG00000212312", "ENSG00000252877", "ENSG00000252623", 
"ENSG00000276036", "ENSG00000200719", "ENSG00000223113", "ENSG00000252029", "ENSG00000212138", "ENSG00000199572", "ENSG00000201413", "ENSG00000222207", "ENSG00000199845", 
"ENSG00000200890", "ENSG00000252951", "ENSG00000201059", "ENSG00000201856", "ENSG00000252999", "ENSG00000212333", "ENSG00000222416", "ENSG00000200028", "ENSG00000276871", 
"ENSG00000251990", "ENSG00000199508", "ENSG00000253080", "ENSG00000202472", "ENSG00000253055", "ENSG00000239021", "ENSG00000272351", "ENSG00000223203", "ENSG00000199592", 
"ENSG00000212276", "ENSG00000200246", "ENSG00000202290", "ENSG00000222832", "ENSG00000252833", "ENSG00000283568", "ENSG00000284736", "ENSG00000252272", "ENSG00000222054", 
"ENSG00000252268", "ENSG00000200343", "ENSG00000212571", "ENSG00000252094", "ENSG00000200114", "ENSG00000212549", "ENSG00000212536", "ENSG00000202187", "ENSG00000201347", 
"ENSG00000223003", "ENSG00000199953", "ENSG00000201966", "ENSG00000222178", "ENSG00000200613", "ENSG00000251890", "ENSG00000253093", "ENSG00000202502", "ENSG00000201923", 
"ENSG00000200370", "ENSG00000252261", "ENSG00000252673", "ENSG00000212237", "ENSG00000222232", "ENSG00000252259", "ENSG00000252970", "ENSG00000212576", "ENSG00000251823", 
"ENSG00000251756", "ENSG00000199337", "ENSG00000252780", "ENSG00000202383", "ENSG00000202047", "ENSG00000201210", "ENSG00000202060", "ENSG00000202422", "ENSG00000275877", 
"ENSG00000200238", "ENSG00000207341", "ENSG00000252598", "ENSG00000201527", "ENSG00000275305", "ENSG00000275215", "ENSG00000274759", "ENSG00000200473", "ENSG00000199806", 
"ENSG00000252479", "ENSG00000201041", "ENSG00000201168", "ENSG00000199994", "ENSG00000222123", "ENSG00000201713", "ENSG00000199809", "ENSG00000252941", "ENSG00000201962", 
"ENSG00000223177", "ENSG00000202164", "ENSG00000200293", "ENSG00000201727", "ENSG00000199638", "ENSG00000199395", "ENSG00000252336", "ENSG00000200558", "ENSG00000251746", 
"ENSG00000222520", "ENSG00000199202", "ENSG00000251785", "ENSG00000252642", "ENSG00000277488", "ENSG00000200839", "ENSG00000251997", "ENSG00000222244", "ENSG00000252086", 
"ENSG00000252289", "ENSG00000200411", "ENSG00000212365", "ENSG00000199900", "ENSG00000252649", "ENSG00000277739", "ENSG00000212542", "ENSG00000222236", "ENSG00000274164", 
"ENSG00000222208", "ENSG00000251760", "ENSG00000200327", "ENSG00000199455", "ENSG00000200313", "ENSG00000200852", "ENSG00000222740", "ENSG00000251857", "ENSG00000252366", 
"ENSG00000199350", "ENSG00000200914", "ENSG00000199620", "ENSG00000251924", "ENSG00000201999", "ENSG00000201523", "ENSG00000201594", "ENSG00000222778", "ENSG00000222182", 
"ENSG00000212628", "ENSG00000200873", "ENSG00000201708", "ENSG00000201274", "ENSG00000200738", "ENSG00000199334", "ENSG00000252959", "ENSG00000199733", "ENSG00000252262", 
"ENSG00000200985", "ENSG00000199318", "ENSG00000252680", "ENSG00000251941", "ENSG00000202054", "ENSG00000199322", "ENSG00000200516", "ENSG00000201492", "ENSG00000272253", 
"ENSG00000274663", "ENSG00000201620", "ENSG00000199874", "ENSG00000251816", "ENSG00000212280", "ENSG00000200619", "ENSG00000222675", "ENSG00000222427", "ENSG00000202356", 
"ENSG00000212171", "ENSG00000252451", "ENSG00000238765", "ENSG00000202526", "ENSG00000252650", "ENSG00000222952", "ENSG00000251993", "ENSG00000251781", "ENSG00000222983", 
"ENSG00000199396", "ENSG00000253083", "ENSG00000212265", "ENSG00000200310", "ENSG00000223013", "ENSG00000251717", "ENSG00000201014", "ENSG00000278233", "ENSG00000200275", 
"ENSG00000253040", "ENSG00000276700", "ENSG00000222578", "ENSG00000222920", "ENSG00000199240", "ENSG00000252637", "ENSG00000252287", "ENSG00000274228", "ENSG00000212527", 
"ENSG00000199354", "ENSG00000252647", "ENSG00000252830", "ENSG00000200079", "ENSG00000201790", "ENSG00000251920", "ENSG00000252618", "ENSG00000201588", "ENSG00000200408", 
"ENSG00000252424", "ENSG00000201321", "ENSG00000251936", "ENSG00000201148", "ENSG00000223262", "ENSG00000201704", "ENSG00000252747", "ENSG00000199839", "ENSG00000199270", 
"ENSG00000274408", "ENSG00000202264", "ENSG00000223290", "ENSG00000199415", "ENSG00000202248", "ENSG00000252041", "ENSG00000276442", "ENSG00000222682", "ENSG00000252866", 
"ENSG00000201096", "ENSG00000252587", "ENSG00000212499", "ENSG00000202063", "ENSG00000201931", "ENSG00000200278", "ENSG00000223086", "ENSG00000252437", "ENSG00000201476", 
"ENSG00000200624", "ENSG00000252267", "ENSG00000252956", "ENSG00000202263", "ENSG00000252070", "ENSG00000200786", "ENSG00000199364", "ENSG00000200381", "ENSG00000200021", 
"ENSG00000202324", "ENSG00000199523", "ENSG00000222383", "ENSG00000200058", "ENSG00000252161", "ENSG00000212505", "ENSG00000199480", "ENSG00000212601", "ENSG00000252535", 
"ENSG00000252828", "ENSG00000201469", "ENSG00000212331", "ENSG00000222418", "ENSG00000199509", "ENSG00000222459", "ENSG00000201285", "ENSG00000252315", "ENSG00000251768", 
"ENSG00000199910", "ENSG00000277411", "ENSG00000222317", "ENSG00000212154", "ENSG00000200036", "ENSG00000199402", "ENSG00000252301", "ENSG00000201325", "ENSG00000252512", 
"ENSG00000277418", "ENSG00000223238", "ENSG00000200468", "ENSG00000200687", "ENSG00000252624", "ENSG00000252553", "ENSG00000199373", "ENSG00000199454", "ENSG00000212396", 
"ENSG00000223076", "ENSG00000278457", "ENSG00000200336", "ENSG00000201035", "ENSG00000222741", "ENSG00000251850", "ENSG00000251729", "ENSG00000222118", "ENSG00000199929", 
"ENSG00000252696", "ENSG00000223019", "ENSG00000201595", "ENSG00000199786", "ENSG00000201861", "ENSG00000274097", "ENSG00000212454", "ENSG00000199985", "ENSG00000223293", 
"ENSG00000202331", "ENSG00000202257", "ENSG00000201185", "ENSG00000251978", "ENSG00000251766", "ENSG00000202147", "ENSG00000223162", "ENSG00000199564", "ENSG00000201876", 
"ENSG00000201728", "ENSG00000252169", "ENSG00000202175", "ENSG00000201415", "ENSG00000252674", "ENSG00000252012", "ENSG00000202322", "ENSG00000275757", "ENSG00000201361", 
"ENSG00000223259", "ENSG00000201763", "ENSG00000202521", "ENSG00000252509", "ENSG00000222205", "ENSG00000202281", "ENSG00000199843", "ENSG00000201695", "ENSG00000252936", 
"ENSG00000201822", "ENSG00000202092", "ENSG00000278189", "ENSG00000239184", "ENSG00000222835", "ENSG00000252342", "ENSG00000201610", "ENSG00000252211", "ENSG00000252516", 
"ENSG00000252902", "ENSG00000251873", "ENSG00000201447", "ENSG00000202225", "ENSG00000278294", "ENSG00000251764", "ENSG00000222921", "ENSG00000222849", "ENSG00000252856", 
"ENSG00000212204", "ENSG00000200487", "ENSG00000199450", "ENSG00000200806", "ENSG00000200587", "ENSG00000201355", "ENSG00000280646", "ENSG00000199407", "ENSG00000222312", 
"ENSG00000201812", "ENSG00000202193", "ENSG00000253031", "ENSG00000222552", "ENSG00000238627", "ENSG00000202386", "ENSG00000252207", "ENSG00000201567", "ENSG00000252164", 
"ENSG00000252376", "ENSG00000201145", "ENSG00000207186", "ENSG00000222302", "ENSG00000252546", "ENSG00000222108", "ENSG00000200601", "ENSG00000199315", "ENSG00000201420", 
"ENSG00000222251", "ENSG00000199352", "ENSG00000283291", "ENSG00000201942", "ENSG00000283433", "ENSG00000252001", "ENSG00000200434", "ENSG00000199804", "ENSG00000252726", 
"ENSG00000251698", "ENSG00000252456", "ENSG00000222346", "ENSG00000222378", "ENSG00000252060", "ENSG00000202334", "ENSG00000283274", "ENSG00000277004", "ENSG00000222585", 
"ENSG00000252714", "ENSG00000212238", "ENSG00000252957", "ENSG00000201518", "ENSG00000252977", "ENSG00000201925", "ENSG00000201939", "ENSG00000202474", "ENSG00000200248", 
"ENSG00000206882", "ENSG00000222428", "ENSG00000252346", "ENSG00000201532", "ENSG00000200227", "ENSG00000199545", "ENSG00000252848", "ENSG00000222608", "ENSG00000212433", 
"ENSG00000238405", "ENSG00000200741", "ENSG00000222838", "ENSG00000252496", "ENSG00000201301", "ENSG00000222285", "ENSG00000252950", "ENSG00000212625", "ENSG00000222419", 
"ENSG00000212595", "ENSG00000212425", "ENSG00000199299", "ENSG00000207277", "ENSG00000274917", "ENSG00000201440", "ENSG00000251705", "ENSG00000199525", "ENSG00000212373", 
"ENSG00000222407", "ENSG00000201000", "ENSG00000222971", "ENSG00000212176", "ENSG00000252905", "ENSG00000200711", "ENSG00000251953", "ENSG00000222955", "ENSG00000238391", 
"ENSG00000222922", "ENSG00000199585", "ENSG00000201312", "ENSG00000274059", "ENSG00000252313", "ENSG00000277049", "ENSG00000212525", "ENSG00000252364", "ENSG00000252246", 
"ENSG00000212251", "ENSG00000251829", "ENSG00000202411", "ENSG00000238908", "ENSG00000201356", "ENSG00000273730", "ENSG00000271924", "ENSG00000222500", "ENSG00000200872"
];
PorcentajeRibo = {};

fout = open("NormalizedExpression.txt", "w");
fout.write("\tID\tName");
for n in NOMBRES:
	PorcentajeMT[n] = 0;
	PorcentajeRibo[n] = 0;
	fout.write("\t"+n);
fout.write("\n");
for i in range(0, len(TABLA["Nombres"])):
	if(CON_NOMBRE_NORMAL): nombreHumano = nombreGen[TABLA["Nombres"][i]];
	else: nombreHumano="";
	fout.write(str(i+1)+"\t"+TABLA["Nombres"][i]+"\t"+nombreHumano);
	for n in NOMBRES:
		fout.write("\t"+str(TABLA[n][i]));
		if TABLA["Nombres"][i] in MT_GENES: PorcentajeMT[n] += TABLA[n][i]
		if TABLA["Nombres"][i] in RIBO_GENES: PorcentajeRibo[n] += TABLA[n][i]
	fout.write("\n");
fout.close();

TABLALOG = {};
TABLALOG["Nombres"] = TABLA["Nombres"];
flog = open("NormalizedExpressionLog.txt", "w");
flog.write("\tID\tName");
for n in NOMBRES:
	flog.write("\t"+n);
	TABLALOG[n] = [];
flog.write("\n");
for i in range(0, len(TABLA["Nombres"])):
	if(CON_NOMBRE_NORMAL): nombreHumano = nombreGen[TABLA["Nombres"][i]];
	else: nombreHumano="";
	flog.write(str(i+1)+"\t"+TABLA["Nombres"][i]+"\t"+nombreHumano);
	for n in NOMBRES:
		aux = log(TABLA[n][i]+1, 2);
		flog.write("\t"+str(aux));
		TABLALOG[n].append(aux);
	flog.write("\n");
flog.close();




fout = open("PercentageMT.txt", "w");
fout.write("Name\tPercentageMT\n");
for n in NOMBRES:
	PorcentajeMT[n] = PorcentajeMT[n]*100.0/sum(TABLA[n]);
	fout.write(n+"\t"+str(PorcentajeMT[n]) + "\n");
fout.close();

fout = open("PercentageRibo.txt", "w");
fout.write("Name\tPercentageRibo\n");
for n in NOMBRES:
	PorcentajeRibo[n] = PorcentajeRibo[n]*100.0/sum(TABLA[n]);
	fout.write(n+"\t"+str(PorcentajeRibo[n]) + "\n");
fout.close();


#Log
for n in NOMBRES:
	for i in range(0, len(TABLA[n])): TABLA[n][i] = log(TABLA[n][i]+1, 2);

################Calculamos columnas especiales
print "AGGREGATE COLUMNS"

Especiales = {
"Time_1": ["CH19075_1", "CH19083_1", "CH19088_1", "CH19091_1", "CH19101_1", "DH16044_1", "DH16076_1", "CH19078_1"],
"Time_2": ["CH16013_2", "CH19083_2", "CH19088_2", "CH19091_2", "CH19101_2", "DH16044_2", "DH16076_2", "CH16023_2"],
"Time_3": ["CH16013_3", "CH19075_3", "CH19083_3", "CH19088_3", "CH19091_3", "CH19101_3", "DH16044_3", "CH19078_3", "CH16023_3"],
"L8": ["DH16044_3", "CH19091_1", "CH19091_2", "CH19091_3", "CH19083_3", "CH19088_1", "CH16023_2", "CH19101_3", "DH16044_1", "CH19075_1", "CH19075_3"],
"L1": ["CH19083_1", "CH19088_2", "CH19078_3", "DH16076_1", "CH16013_2", "CH16023_3","CH19101_1", "DH16044_2", "CH19083_2", "CH19088_3", "CH19078_1", "DH16076_2", "CH16013_3", "CH19101_2"],
}


ELEMENTOS_DE = {}
ELEMENTOS_DE_LOG = {}

numExluidos = 0;
for e in Especiales:
	TABLA[e] = [];
	ELEMENTOS_DE[e] = {};
	ELEMENTOS_DE_LOG[e] = {};
	for i in range(0, len(TABLA["Nombres"])):
		aux = [];
		aux2 = [];
		for n in Especiales[e]:
			aux.append((2**TABLA[n][i])-1);
			aux2.append(TABLA[n][i]);
				
		TABLA[e].append(log(np.mean(aux)+1, 2));
		ELEMENTOS_DE[e][TABLA["Nombres"][i]] = aux;
		ELEMENTOS_DE_LOG[e][TABLA["Nombres"][i]] = aux2;

print numExluidos

##############\ Calculamos columnas especiales
print "COMPARATIVES"

COMPARATIVAS = [];
PEARSONS = [];
SPEARMANS = [];
for n in nombres:
	auxP = []
	auxS = [];
	for n2 in nombres:
		aux = [n, n2]
		COMPARATIVAS.append(aux);
		auxP.append(pearsonr(TABLA[n], TABLA[n2])[0]);
		auxS.append(spearmanr(TABLA[n], TABLA[n2])[0]);
	PEARSONS.append(auxP);
	SPEARMANS.append(auxS);

plt.figure();
plt.pcolor(PEARSONS);
plt.title("Pearson Correlation");
plt.xticks(np.arange(0.5, len(nombres), 1.0),nombres, rotation=27);
plt.yticks(np.arange(0.5, len(nombres), 1.0),nombres);
plt.gca().invert_yaxis();
plt.colorbar();
plt.savefig("Figures/PearsonR.png");


plt.figure();
plt.pcolor(SPEARMANS);
plt.title("Spearman Correlation");
plt.xticks(np.arange(0.5, len(nombres), 1.0),nombres, rotation=27);
plt.yticks(np.arange(0.5, len(nombres), 1.0),nombres);
plt.gca().invert_yaxis();
plt.colorbar();
plt.savefig("Figures/SpearmanR.png");
#~ plt.show();
#~ exit();
k = 0;

COMPARATIVAS = [
["Time_1", "Time_2"],
["Time_1", "Time_3"],
["Time_3", "Time_2"],
["L8", "L1"]
]

for c in COMPARATIVAS:
	
	plt.clf();
	
	plt.title(c[0] + " vs " + c[1]);
	plt.scatter(TABLA[c[0]], TABLA[c[1]],color="blue", alpha=0.7);
	plt.xlabel("Log2("+c[0]+")");
	plt.ylabel("Log2("+c[1]+")");
	plt.xlim(-0.1, 20);
	plt.ylim(-0.1, 20);
	plt.savefig("Figures/"+c[0]+"_VS_"+c[1]+".png");

	


################## DEG Analysis

FCLimit = log(1.5, 2);
pValueLimit = 0.05

for c in COMPARATIVAS:

	print c

	pValues = [];
	FCs = [];
	pValuesRepresentar = [];

	if CON_NOMBRE_NORMAL: ENCABEZADO = "Transcipt ID\tGene name\tLog2("+c[0]+")\tLog2("+c[1]+")\tLog2(FC)\tpValue\tFDR < 0.05\n"
	else: ENCABEZADO = "Transcipt ID\tLog2("+c[0]+")\tLog2("+c[1]+")\tLog2(FC)\tpValue\tFDR < 0.05\n"

	for t in transcritos:
		Sano = log(np.mean(ELEMENTOS_DE[c[0]][t])+1,2);
		Enf = log(np.mean(ELEMENTOS_DE[c[1]][t])+1,2);
		DE = Enf-Sano;
		pValue = ttest_ind(ELEMENTOS_DE[c[0]][t], ELEMENTOS_DE[c[1]][t])[1]
		if np.isnan(pValue): pValue=1.0
		
		pValuesRepresentar.append(-log(pValue, 10));
		pValues.append(pValue);
		FCs.append(DE);
		
	fout.close();

	plt.clf()
	plt.title("Uncorrected vPlot");
	plt.scatter(FCs, pValuesRepresentar, alpha=0.5);
	plt.xlabel("Log2(FC)");
	plt.ylabel("-Log10(pValue)");
	plt.xlim(-6, 6);
	plt.ylim(0, 10);
	plt.savefig("Figures/UncorrectedVolcano_"+c[0]+"_VS_"+c[1]+".png");


	fUp = open("Overexpressed_"+c[0]+"_VS_"+c[1]+".txt", "w");
	fDown = open("Repressed_"+c[0]+"_VS_"+c[1]+".txt", "w");
	fout = open(c[0]+"_VS_"+c[1]+"_All_Corrected.txt", "w");
	fout.write(ENCABEZADO);
	fUp.write(ENCABEZADO);
	fDown.write(ENCABEZADO);

	colors=[]
	pValuesRepresentar=[];
	i = 0;
	CorrectedpValue = multipletests(pValues, 0.05, method='fdr_bh')[1];
	for t in transcritos:
		Sano = log(np.mean(ELEMENTOS_DE[c[0]][t])+1,2);
		Enf = log(np.mean(ELEMENTOS_DE[c[1]][t])+1,2);
		pValue=pValues[i];
		pValueCorregido=CorrectedpValue[i]
		FC = FCs[i]
		pValuesRepresentar.append(-log(CorrectedpValue[i],10))
		colors.append("blue");
		if CON_NOMBRE_NORMAL:linea = t+"\t"+nombreGen[t]+"\t"+str(Sano)+"\t"+str(Enf)+"\t"+str(FC)+"\t"+str(pValue)+"\t"+str(pValueCorregido)+"\n";
		else:linea = t+"\t"+str(Sano)+"\t"+str(Enf)+"\t"+str(FC)+"\t"+str(pValue)+"\t"+str(pValueCorregido)+"\n";
		if(pValueCorregido <= pValueLimit):
			H = np.mean(ELEMENTOS_DE[c[0]][t]);
			P = np.mean(ELEMENTOS_DE[c[1]][t]);
			if(FC > FCLimit): #UP
				if(H >= LIMITE_EXPR and P >= LIMITE_EXPR):
					fUp.write(linea);
					colors[i] = "red";
				elif(P >= 2*LIMITE_EXPR):
					fUp.write(linea);
					colors[i] = "red";
			elif(FC < -FCLimit): #Down
				if(H >= LIMITE_EXPR and P >= LIMITE_EXPR):
					fDown.write(linea);
					colors[i] = "red";
				elif(H >= 2*LIMITE_EXPR):
					fDown.write(linea);
					colors[i] = "red";
		fout.write(linea);
		i+=1;
	fout.close();
	fUp.close();
	fDown.close();

	plt.clf()
	plt.title("vPlot");
	plt.scatter(FCs, pValuesRepresentar, color=colors, alpha=0.5);
	plt.xlabel("Log2(FC)");
	plt.ylabel("-Log10(pValue)");
	plt.xlim(-6, 6);
	plt.ylim(0, 5);
	plt.savefig("Figures/CorrectedVolcano_"+c[0]+"_VS_"+c[1]+".png");

	plt.clf();
	plt.hist(pValues, np.arange(0,1,0.01));
	plt.savefig("Figures/Uncorrected_pValues_"+c[0]+"_VS_"+c[1]+".png");


	plt.clf();
	plt.hist(CorrectedpValue, np.arange(0,1,0.01));
	plt.savefig("Figures/Corrected_pValues_"+c[0]+"_VS_"+c[1]+".png");

print "DONE"
