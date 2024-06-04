#Close terminal
ssh serv@000.000.00.000
password: XXXXXXXXXXXXX

#Downloads (NCBI)
fastq-dump SRR12899120 -v --gzip
fastq-dump SRR12899121 -v --gzip
fastq-dump SRR12899122 -v --gzip
fastq-dump SRR12899123 -v --gzip
fastq-dump SRR27004683 -v --gzip
fastq-dump SRR27004682 -v --gzip
fastq-dump SRR27004681 -v --gzip
fastq-dump SRR23881846 -v --gzip
fastq-dump SRR23881845 -v --gzip
fastq-dump SRR23881844 -v --gzip
fastq-dump SRR11906443 -v --gzip
fastq-dump SRR11906444 -v --gzip
fastq-dump SRR11906445 -v --gzip
fastq-dump SRR11906446 -v --gzip
fastq-dump SRR11906447 -v --gzip
fastq-dump SRR11906448 -v --gzip
fastq-dump SRR15602503 -v --gzip
fastq-dump SRR15602502 -v --gzip
fastq-dump SRR15602501 -v --gzip
fastq-dump SRR15602500 -v --gzip
fastq-dump SRR12899124 -v --gzip
fastq-dump SRR12899125 -v --gzip
fastq-dump SRR12899126 -v --gzip
fastq-dump SRR12899127 -v --gzip
fastq-dump SRR27004680 -v
fastq-dump SRR27004679 -v
fastq-dump SRR27004678 -v
fastq-dump SRR23881852 -v
fastq-dump SRR23881851 -v
fastq-dump SRR23881850 -v
fastq-dump SRR11906450 -v
fastq-dump SRR11906451 -v
fastq-dump SRR11906452 -v
fastq-dump SRR11906453 -v
fastq-dump SRR11906454 -v
fastq-dump SRR11906455 -v
fastq-dump SRR15602499 -v
fastq-dump SRR15602498 -v
fastq-dump SRR15602497 -v
fastq-dump SRR15602496 -v
fastq-dump SRR11657441 -v
fastq-dump SRR11657442 -v
fastq-dump SRR11657443 -v
fastq-dump SRR11657429 -v
fastq-dump SRR11657430 -v
fastq-dump SRR11657431 -v
fastq-dump SRR776573 -v
fastq-dump SRR776577 -v
fastq-dump SRR776569 -v
fastq-dump SRR776557 -v
fastq-dump SRR776561 -v
fastq-dump SRR776553 -v
fastq-dump SRR12145330 -v
fastq-dump SRR12145331 -v
fastq-dump SRR12145332 -v
fastq-dump SRR12145334 -v
fastq-dump SRR12145337 -v
fastq-dump SRR12145338 -v
fastq-dump SRR21987044 -v
fastq-dump SRR21987045 -v
fastq-dump SRR21987046 -v
fastq-dump SRR21987038 -v
fastq-dump SRR21987039 -v
fastq-dump SRR21987040 -v
fastq-dump SRR4228541 -v
fastq-dump SRR4228542 -v
fastq-dump SRR4228543 -v
fastq-dump SRR4228545 -v
fastq-dump SRR4228546 -v
fastq-dump SRR4228547 -v
fastq-dump SRR18713679 -v
fastq-dump SRR18713680 -v
fastq-dump SRR18713681 -v
fastq-dump SRR18713675 -v
fastq-dump SRR18713676 -v
fastq-dump SRR18713678 -v
fastq-dump SRR17646267 -v
fastq-dump SRR17646268 -v
fastq-dump SRR17646269 -v
fastq-dump SRR17646264 -v
fastq-dump SRR17646265 -v
fastq-dump SRR17646266 -v
fastq-dump SRR16538453 -v
fastq-dump SRR16538454 -v
fastq-dump SRR16538455 -v
fastq-dump SRR16538444 -v
fastq-dump SRR16538445 -v
fastq-dump SRR16538446 -v
fastq-dump SRR8309415 -v
fastq-dump SRR8309416 -v
fastq-dump SRR8309417 -v
fastq-dump SRR8309418 -v
fastq-dump SRR8309419 -v
fastq-dump SRR8309420 -v
ls






#mv File path 

#Create directories for files
mkdir Canis  Homo  Mus  Oryctolagus  Vespertilio #para comportar as amostras separadamente por organismo
mv "/Downloads/files" /file/directories/
mkdir raw_reads #Para comportar os diretórios com dados brutos
mv Canis  Homo  Mus  Oryctolagus  Vespertilio /media/koch/hd1_4tb/samuelR/Gil_Project/raw_reads #para agrupar os diretórios

#Create output directories
cd /media/koch/hd1_4tb/samuelR/Gil_Project/ #se não estiver no diretório Gil_Project
mkdir data #para comportar os outputs
mkdir tools #caso precise organizar os programas em um local conhecido

mkdir fastqc_output prinseq_output salmon_output
mv fastqc_output prinseq_output salmon_output /media/koch/hd1_4tb/samuelR/Gil_Project/data

#Checking if the tools are installed
fastqc --version #ok FastQC v0.11.9
prinseq-lite --version #ok PRINSEQ-lite 0.20.4
salmon --version #ok salmon 1.4.0

#Data quality analysis
#Run data - FASTQC
fastqc arquivo.fastq -t 3 #Realizar o carregamento individualmente para cada arquivo
rm files_fastqc.zip #Opcional
mv files_fastqc.html /media/koch/hd1_4tb/samuelR/Gil_Project/data/fastqc_output
/media/koch/hd1_4tb/samuelR/Gil_Project/data/fastqc_output mkdir Canis_fastqc  Homo_fastqc  Mus_fastqc  Oryctolagus_fastqc  Vespertilio_fastqc
/media/koch/hd1_4tb/samuelR/Gil_Project/data/fastqc_output mv files_fastqc.html /media/koch/hd1_4tb/samuelR/Gil_Project/data/fastqc_output/Organismo_fastqc

#View HTML files (close a new terminal)
mkdir HTMLs #para receber os arquivos copiados do servidor na máquina
scp koch@200.239.83.194:/media/koch/hd1_4tb/samuelR/Gil_Project/data/fastqc_output/Canis_fastqc/FILE_fastqc.html . #para copiar os arquivos do serv

#Run data - PRINSEQ (raw_reads)
prinseq-lite -fastq FILE.fastq -min_qual_mean 20 -out_good FILE_filtered_good -out_bad null 
mv FILE_filtered_good /media/koch/hd1_4tb/samuelR/Gil_Project/data/prinseq_output/Organismo_prinseq #Dentro de prinseq_output é indicado criar novos diretórios (por organismo) para melhor organizar os outputs



#Run data - Salmon (antes de fazer a predição é necessária a criação de um 'index' com um genoma de referência que pode ser baixado do Ensembl em formato .fasta)
salmon index -t Homo_sapiens.fa -i Homo_index #para criação do index que gera o diretório Homo_index com todos os arquivos necessários para a próx etapa
salmon quant -i caminho_do_index_criado -l A -r /caminho/da/triplicata/01.sf -o salmon_output


-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#Todos os arquivos de output do salmon serão nomeados automaticamente como "quant.sf" e podem ser renomeados com os nomes das sequencias e agrupados em um diretório só. Ex (SRR111.sf, SRR222.sf, SRR333.sf).
#É interessante que eles sejam renomeados para que as duas triplicatas (controle e condição) possam ser agrupadas no mesmo diretório, facilitando o carregamento para o ambiente R.
#Os passos seguintes acontecerão em ambiente Rstudio.



#Run data - IsoformSwitchAnalyzer (Bioconductor, R).
#Importing the base directory

# Load the IsoformSwitchAnalyzeR library
library(IsoformSwitchAnalyzeR)

# Definir diretórios e arquivos
quantsDir <- "/home/samuel/Desktop/DISCIPLINAS/GIL_02/PROJECT/WF/quant_sf/Homo"
gtfDir <- "/home/samuel/Desktop/DISCIPLINAS/GIL_02/PROJECT/WF/IMPORTANTES/Homo_sapiens.GRCh38.112.chr.gtf"
cdnaDir <- "/home/samuel/Desktop/DISCIPLINAS/GIL_02/PROJECT/WF/IMPORTANTES/Homo_sapiens.GRCh38.cdna.all.fa"
resDir <- "/home/samuel/Desktop/DISCIPLINAS/GIL_02/PROJECT/WF/final"

# Definir design matrix
design.df <- data.frame(
  sampleID = c("Sample_1", "Sample_2", "Sample_3", "Sample_4", "Sample_5", "Sample_6"), 
  condition = c("Control_1", "Control_1", "Control_1", "Condition_1", "Condition_1", "Condition_1")
)

# Definir comparação
comp.df <- data.frame(
  condition_1 = c("Control_1"), 
  condition_2 = c("Condition_1")
)

# Vetor de amostras com caminhos completos
sampleVector <- c(
  "/home/samuel/Desktop/DISCIPLINAS/GIL_02/PROJECT/WF/quant_sf/Homo/Sample_1/SRR23881844.sf",
  "/home/samuel/Desktop/DISCIPLINAS/GIL_02/PROJECT/WF/quant_sf/Homo/Sample_2/SRR23881845.sf",
  "/home/samuel/Desktop/DISCIPLINAS/GIL_02/PROJECT/WF/quant_sf/Homo/Sample_3/SRR23881846.sf",
  "/home/samuel/Desktop/DISCIPLINAS/GIL_02/PROJECT/WF/quant_sf/Homo/Sample_4/SRR23881850.sf",
  "/home/samuel/Desktop/DISCIPLINAS/GIL_02/PROJECT/WF/quant_sf/Homo/Sample_5/SRR23881851.sf",
  "/home/samuel/Desktop/DISCIPLINAS/GIL_02/PROJECT/WF/quant_sf/Homo/Sample_6/SRR23881852.sf"
)

# Importar dados de expressão de isoformas
quants.df <- importIsoformExpression(
  sampleVector = sampleVector,
  addIsofomIdAsColumn = TRUE
)

# Ajustar nomes das colunas de counts e abundance para corresponder aos IDs de amostras, exceto a primeira coluna
colnames(quants.df$counts)[-1] <- design.df$sampleID
colnames(quants.df$abundance)[-1] <- design.df$sampleID

# Verificar a estrutura após renomeação
print(head(quants.df$counts))
print(head(quants.df$abundance))

# Importar dados para análise
aSwitchlist <- importRdata(
  isoformCountMatrix = quants.df$counts, 
  isoformRepExpression = quants.df$abundance, 
  designMatrix = design.df, 
  isoformExonAnnoation = gtfDir, 
  isoformNtFasta = cdnaDir, 
  comparisonsToMake = comp.df,
  fixStringTieAnnotationProblem = TRUE
)

# Pré-filtragem dos dados
aSwitchlist <- preFilter(
  aSwitchlist,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  IFcutoff = 0,
  removeSingleIsoformGenes = TRUE,
  reduceToSwitchingGenes = FALSE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  onlySigIsoforms = TRUE,
  keepIsoformInAllConditions = FALSE,
  quiet = FALSE
)

# Teste de comutação de isoformas usando DEXSeq
aSwitchlist <- isoformSwitchTestDEXSeq(
  aSwitchlist, 
  reduceToSwitchingGenes = FALSE
)

# Análise combinada de comutação de isoformas
aSwitchlist <- isoformSwitchAnalysisCombined(
  aSwitchlist,
  pathToOutput = resDir 
)

# Análise de splicing alternativo
aSwitchlist <- analyzeAlternativeSplicing(
  aSwitchlist,
  onlySwitchingGenes = FALSE
)

# Análise de retenção de íntrons
aSwitchlist <- analyzeIntronRetention(
  aSwitchlist, 
  onlySwitchingGenes = FALSE
)

# Função para extrair e salvar gráficos e tabelas
extractPlots <- function(){
  pdf(paste0(resDir, "/outputPlots.pdf"))
  extractSplicingSummary(aSwitchlist, plotGenes = TRUE, returnResult = FALSE)
  extractSplicingSummary(aSwitchlist, plotGenes = FALSE, returnResult = FALSE)
  extractSplicingEnrichment(aSwitchlist, returnResult = FALSE)
  extractConsequenceSummary(aSwitchlist)
  extractSwitchOverlap(aSwitchlist, plotIsoforms = TRUE, plotSwitches = FALSE, plotGenes = FALSE)
  extractSwitchOverlap(aSwitchlist, plotIsoforms = FALSE, plotSwitches = TRUE, plotGenes = FALSE)
  extractSwitchOverlap(aSwitchlist, plotIsoforms = FALSE, plotSwitches = FALSE, plotGenes = TRUE)
  dev.off()
  
  write.csv(extractSplicingSummary(aSwitchlist, plot = FALSE, returnResult = TRUE), file = paste0(resDir, "/splicingSummary.csv"))
  write.csv(extractSplicingEnrichment(aSwitchlist, plot = FALSE, returnResult = TRUE), file = paste0(resDir, "/splicingEnrichment.csv"))
  write.csv(extractSwitchSummary(aSwitchlist), file = paste0(resDir, "/switchSummary.csv"))
  write.csv(extractTopSwitches(aSwitchlist, n = NA, inEachComparison = TRUE), file = paste0(resDir, "/topSwitches.csv"))
}

# Executar a função para extrair e salvar resultados
extractPlots()