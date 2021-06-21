def lectura_porcr (cromosoma, archivo,  cromosomas,posiciones_crom,tipo_archivo): #Va leyendo cromosoma a cromosoma, empezando por donde lo dejo el anterior
	with open(archivo, "r") as archivo:
		archivo_read = []
		archivo.seek(posiciones_crom[tipo_archivo][cromosoma]) #Busca la posicion del cromosoma que vamos a leer
		line = archivo.readline()
		posicion_linea = posiciones_crom[tipo_archivo][cromosomas[cromosomas.index(cromosoma)]]
		bool_linea = 0
		bool_cond = 0
		while line:
			if line[0:len(cromosoma)+1] != (cromosoma + "\t"): #Comprueba si hemos cambiado de cromosoma, si es así, guarda la posicion
				bool_cond = 1
				try:
					posiciones_crom[tipo_archivo][cromosomas[cromosomas.index(cromosoma)+1]] = posicion_linea
				except IndexError:
					pass
				break
			else:
				archivo_read.append(line)
				posicion_linea = archivo.tell() #Guarda la posición de la última línea
				bool_linea = 1
				line = archivo.readline()
	if bool_cond == 0: #Si sale sin pasar por el if, se ha acabado el archivo
		try:
			posiciones_crom[tipo_archivo][cromosomas[cromosomas.index(cromosoma)+1]] = posiciones_crom[tipo_archivo][cromosomas[cromosomas.index(cromosoma)]]
		except IndexError:
			pass
	if bool_linea == 0: #Si sale sin pasar por las lineas, devuelve nada
		try:
			posiciones_crom[tipo_archivo][cromosomas[cromosomas.index(cromosoma)+1]] = posicion_linea
		except IndexError:
			pass
		return []
	else:
		return archivo_read #Devuelve las lineas sin separarlas  
def check_gene (align_set, inicio_seq_alin, final_seq_alin, inicio_seq_gen, final_seq_gen, tantox1_identidad):
    min_longitud = tantox1_identidad*(len(align_set)) 
    #Compruebo condiciones para aumentar el rendimiento
    if final_seq_alin < inicio_seq_gen: #La anotacion esta por delante del alineamiento, hay que salir de las anotaciones y avanzar los alineamientos
        return 0
    elif inicio_seq_alin > final_seq_gen: #El alineamiento esta por delante de la anotacion, se pasa al siguiente alineamiento
        return 2
    if inicio_seq_alin >= inicio_seq_gen: #Compruebo cada posible condicion para que el alineamiento caiga en la anotacion
        inicio_overlap = inicio_seq_alin
    else:
        inicio_overlap = inicio_seq_gen
    if final_seq_alin <= final_seq_gen:
        final_overlap = final_seq_alin
    else:
        final_overlap = final_seq_gen
    overlap_set = set(tuple(range(inicio_overlap,final_overlap+1)))
    if len(overlap_set) >= min_longitud: #Compruebo que pase del % de longitud deseado
        return 1
def check_hierarchy (tipo, tipo_align):
    hierarchy=["tRNA","Mt_tRNA","snRNA","sncRNA","snoRNA",
    "miRNA","pseudogene","polymorphic_pseudogene","rmsk" ,"protein", "protein_coding","lincRNA","lncRNA", "antisense","sense_intronic",
    "sense_overlapping","3prime_overlapping_ncrna","misc_RNA","processed_transcript"]
    if hierarchy.index(tipo) < hierarchy.index(tipo_align):
        return tipo
    else:
        return tipo_align
def assign_genes(annot_file,align_file,cromosomas,posiciones_crom,strands,low_hierarchy,read_counts,tantox1_identidad,basename):
    for cromosoma in cromosomas:
        print (cromosoma)
        gene_anot = [gene.split() for gene in lectura_porcr(cromosoma, annot_file,cromosomas,posiciones_crom,"annot")]
        genes_fw = [gene_fw for gene_fw in gene_anot if gene_fw[5] == "+" ]
        genes_rv = [gene_rv for gene_rv in gene_anot if gene_rv[5] == "-"]
        align = [alignment.split() for alignment in lectura_porcr(cromosoma, align_file,  cromosomas, posiciones_crom, "align")]
        align_fw = [align_fw for align_fw in align if align_fw[5] == "+"]
        align_rv = [align_rv for align_rv in align if align_rv[5] == "-"]
        for strand in strands:
            print(strand)
            if strand=="fw":
                alignments=align_fw
                genes=genes_fw
            else:
                alignments=align_rv
                genes=genes_rv
            for alignment in alignments:
                inicio_seq_alin =int(alignment[1])
                final_seq_alin =int(alignment[2])
                align_set= set(tuple(range(inicio_seq_alin, final_seq_alin+1)))
                reads = int(alignment[4])
                tipo_align = ""
                for gen in genes:
                    inicio_seq_gen = int(gen[1])
                    final_seq_gen = int(gen[2])
                    tipo = gen[3]
                    check = check_gene (align_set, inicio_seq_alin, final_seq_alin, inicio_seq_gen, final_seq_gen, tantox1_identidad)
                    if check == 0:
                        break
                    elif check == 2:
                        continue
                    elif check ==1:
                        if tipo in ["rRNA","Mt_rRNA"]:
                            tipo_align = tipo
                            break
                        if tipo_align != "":
                            if tipo_align != tipo: 
                                tipo_align = check_hierarchy (tipo, tipo_align)
                        else:
                            tipo_align = tipo
                        if tipo_align in ["rRNA","Mt_rRNA"]:
                            break

                if tipo_align == "":
                    read_counts[6]+= reads
                else:
                    if tipo_align in ["rRNA","Mt_rRNA"]:
                        read_counts[0] += reads
                    elif tipo_align in ["tRNA","Mt_tRNA"]:
                        read_counts[1] += reads
                    elif tipo_align in ["snRNA","snoRNA","miRNA", "scRNA"]:
                        read_counts[2] += reads
                    elif tipo_align in ["pseudogene","polymorphic_pseudogene"]:
                        read_counts[3] += reads
                    elif tipo_align in ["protein_coding"]:
                        read_counts[4] += reads
                    elif tipo_align in ["lincRNA","antisense","sense_intronic","sense_overlapping","3prime_overlapping_ncrna","misc_RNA","processed_transcript"]:
                        read_counts[5] += reads

    with open ("Results_"+str(tantox1_identidad)+"_"+ basename + ".txt", "a") as resultados:
        keys = low_hierarchy
        f_line = "\n----------"
        for key in keys:
            f_line += "\t" + key 
        resultados.write(f_line)
        results = "\n" + "Reads"
        for key in keys:
            results += "\t" + str(read_counts[low_hierarchy.index(key)])
        resultados.write(results)
                    
for a in [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.00000000000001]:
    for b in ["sortedRep2.bed","sortedRep1.bed"]:
        if b=="sortedRep2.bed":
            basename = "Rep2"
        else:
            basename ="Rep1"
        annot_file="Genome.bed"
        align_file=b
        tantox1_identidad=a
        cromosomas = ["chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
        "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6",
        "chr7","chr8","chr9","chrM","chrX","chrY"]
        posiciones_crom={"align":{"chr1":0},
                        "annot":{"chr1":0}}
        strands = ["fw", "rv"]
        low_hierarchy = ["rRNA", "tRNA", "sncRNA", "pseudogen","mRNA","lncRNA", "Intergenico"]
        read_counts=[0 for eslabon in low_hierarchy]
        assign_genes(annot_file,align_file,cromosomas,posiciones_crom,strands,low_hierarchy,read_counts,tantox1_identidad,basename)


