library(wlab.block)
library(data.table)
## caculate interchain distance   K13
interchain_distance_K13 <- bio3D_minimum_interchain_distances(
  input_file = "6h46.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_K13)[2:3]<-paste0(colnames(interchain_distance_K13)[2:3],"_K13")
colnames(interchain_distance_K13)[1] <- c("Pos_real")
anno<-interchain_distance_K13
anno[, K13_bind := 'no']
anno[scHAmin_ligand_K13 <5 , K13_bind := 'yes']
anno[K13_bind=="yes",site_type:="Binding interface site"]
names(anno)
library(dplyr)
binding_positions <- anno %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)
print(binding_positions)



library(wlab.block)
library(data.table)
## caculate interchain distance  K19
interchain_distance_K19 <- bio3D_minimum_interchain_distances(
  input_file = "6h47.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_K19)[2:3]<-paste0(colnames(interchain_distance_K19)[2:3],"_K19")
colnames(interchain_distance_K19)[1] <- c("Pos_real")
anno<-interchain_distance_K19
anno[, K19_bind := 'no']
anno[scHAmin_ligand_K19 <5 , K19_bind := 'yes']
anno[K19_bind=="yes",site_type:="Binding interface site"]
names(anno)
library(dplyr)
binding_positions <- anno %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)
print(binding_positions)




library(wlab.block)
library(data.table)
## caculate interchain distance  K27
interchain_distance_K27 <- bio3D_minimum_interchain_distances(
  input_file = "5o2s.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_K27)[2:3]<-paste0(colnames(interchain_distance_K27)[2:3],"_K27")
colnames(interchain_distance_K27)[1] <- c("Pos_real")
anno<-interchain_distance_K27
anno[, K27_bind := 'no']
anno[scHAmin_ligand_K27 <5 , K27_bind := 'yes']
anno[K27_bind=="yes",site_type:="Binding interface site"]
names(anno)
library(dplyr)
binding_positions <- anno %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)
print(binding_positions)







library(wlab.block)
library(data.table)
## caculate interchain distance  K55
interchain_distance_K55 <- bio3D_minimum_interchain_distances(
  input_file = "5mla.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_K55)[2:3]<-paste0(colnames(interchain_distance_K55)[2:3],"_K55")
colnames(interchain_distance_K55)[1] <- c("Pos_real")
anno<-interchain_distance_K55
anno[, K55_bind := 'no']
anno[scHAmin_ligand_K55 <5 , K55_bind := 'yes']
anno[K55_bind=="yes",site_type:="Binding interface site"]
names(anno)
library(dplyr)
binding_positions <- anno %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)
print(binding_positions)




library(wlab.block)
library(data.table)
library(dplyr)
## caculate interchain distance  RAF1
interchain_distance_RAF1 <- bio3D_minimum_interchain_distances(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6vjj.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_RAF1)[2:3]<-paste0(colnames(interchain_distance_RAF1)[2:3],"_RAF1")
colnames(interchain_distance_RAF1)[1] <- c("Pos_real")
anno<-interchain_distance_RAF1
anno[, RAF1_bind := 'no']
anno[scHAmin_ligand_RAF1 <5 , RAF1_bind := 'yes']
anno[RAF1_bind=="yes",site_type:="Binding interface site"]
#names(anno)
binding_positions <- anno %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)
print(binding_positions)


library(wlab.block)
library(data.table)
## caculate interchain distance  PI3KCG
interchain_distance_PI3KCG <- bio3D_minimum_interchain_distances(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1he8.pdb",
  chain_query = "B",
  chain_target = "A")
colnames(interchain_distance_PI3KCG)[2:3]<-paste0(colnames(interchain_distance_PI3KCG)[2:3],"_PI3KCG")
colnames(interchain_distance_PI3KCG)[1] <- c("Pos_real")
anno<-interchain_distance_PI3KCG
anno[, PI3KCG_bind := 'no']
anno[scHAmin_ligand_PI3KCG <5 , PI3KCG_bind := 'yes']
anno[PI3KCG_bind=="yes",site_type:="Binding interface site"]
#names(anno)
library(dplyr)
binding_positions <- anno %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)
print(binding_positions)






library(wlab.block)
library(data.table)
## caculate interchain distance  RALGDS
interchain_distance_RALGDS <- bio3D_minimum_interchain_distances(
  input_file = "1lfd.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_RALGDS)[2:3]<-paste0(colnames(interchain_distance_RALGDS)[2:3],"_RALGDS")
colnames(interchain_distance_RALGDS)[1] <- c("Pos_real")
anno<-interchain_distance_RALGDS
anno[, RALGDS_bind := 'no']
anno[scHAmin_ligand_RALGDS <5 , RALGDS_bind := 'yes']
anno[RALGDS_bind=="yes",site_type:="Binding interface site"]
names(anno)
library(dplyr)
binding_positions <- anno %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)
print(binding_positions)







library(wlab.block)
library(data.table)
## caculate interchain distance  SOS1
interchain_distance_SOS1 <- bio3D_minimum_interchain_distances(
  input_file = "1bkd.pdb",
  chain_query = "R",
  chain_target = "S")
colnames(interchain_distance_SOS1)[2:3]<-paste0(colnames(interchain_distance_SOS1)[2:3],"_SOS1")
colnames(interchain_distance_SOS1)[1] <- c("Pos_real")
anno<-interchain_distance_SOS1
anno[, SOS1_bind := 'no']
anno[scHAmin_ligand_SOS1 <5 , SOS1_bind := 'yes']
anno[SOS1_bind=="yes",site_type:="Binding interface site"]
names(anno)
library(dplyr)
binding_positions <- anno %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)
print(binding_positions)


library(wlab.block)
library(data.table)
## caculate distance from chain to a specific HETATM
chain_hetatm_dis_RAF1 <- bio3D_chain_hetatm_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_RAF1)[2:3]<-paste0(colnames(chain_hetatm_dis_RAF1)[2:3],"_GTP")
colnames(chain_hetatm_dis_RAF1)[1] <- c("Pos_real")
anno<-chain_hetatm_dis_RAF1
anno[, GTP_bind := 'no']
anno[scHAmin_ligand_GTP <5 , GTP_bind := 'yes']
anno[GTP_bind=="yes",site_type:="GTP Binding interface site"]
names(anno)
library(dplyr)
binding_positions <- anno %>%
  filter(site_type == "GTP Binding interface site") %>%
  pull(Pos_real)
print(binding_positions)




library(wlab.block)
library(data.table)
## caculate interchain distance   GDP结合口袋位点
chain_hetatm_dis_K13 <- bio3D_chain_hetatm_distances(
  input_file = "6h46.pdb",
  chain_query = "A",
  hetatm_target = "GDP"
)
colnames(chain_hetatm_dis_K13)[2:3]<-paste0(colnames(chain_hetatm_dis_K13)[2:3],"_GDP")
colnames(chain_hetatm_dis_K13)[1] <- c("Pos_real")
anno<-chain_hetatm_dis_K13
anno[, GDP_bind := 'no']
anno[scHAmin_ligand_GDP <5 , GDP_bind := 'yes']
anno[GDP_bind=="yes",site_type:="Binding interface site"]
names(anno)
library(dplyr)
binding_positions <- anno %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)
print(binding_positions)







