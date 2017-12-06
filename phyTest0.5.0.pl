#!/usr/bin/perl

=cut

Title: phyTest

Version: 0.5.0 (alpha)

Author: Lucas Marques

=cut

use strict;
use warnings;
use List::Util "sum";
use List::MoreUtils "uniq";
#use Math::Gauss; # subrotinas inclusas embaixo para o caso do pacote não estar instalado


# ==================== INICIO DO PROGRAMA ==================== #

#recebe argumentos e identifica, pelos flags, os nomes de arquivos e configurações necessárias aos testes
my $seqs;
my $submodel;
my $treez;
my $nreps = 1000;
my @testconf = ("BP","KH","SH");
my $verb = 0;
my $ncores = 2;
my $redo = "-redo";
print "Reading input...\n";
for(my $i=0; $i<scalar(@ARGV); $i+=2){
	if($ARGV[$i] eq "-s"){
		$seqs = $ARGV[$i+1];
		print "\tSequence alignment file set to '$seqs'\n";
	}
	elsif($ARGV[$i] eq "-m"){
		$submodel = $ARGV[$i+1];
		print "\tSubstitution model set to '$submodel'\n";
	}
	elsif($ARGV[$i] eq "-z"){
		$treez = $ARGV[$i+1];
		print "\tTrees file set to '$treez'\n";
	}
	elsif($ARGV[$i] eq "-n"){
		$nreps = $ARGV[$i+1];
		print "\tNumber of replicates set to '$nreps'\n";
	}
	elsif($ARGV[$i] eq "-t"){
		@testconf = split("/",$ARGV[$i+1]);
		print "\tTests set to '$ARGV[$i+1]'\n";
	}
	elsif($ARGV[$i] eq "-verb"){
		if($ARGV[$i+1] =~ /^T(RUE)?/){
			$verb = 1;
		}
		elsif($ARGV[$i+1] !~ /^F(ALSE)?/){
			die "Invalid entry '$ARGV[$i+1]' for verbose mode (-verb). It must be a boolean (TRUE/FALSE, or shortly T/F). Please, reset and try again.\n";
		}
	}
	elsif($ARGV[$i] eq "-nc"){
		if($ARGV[$i+1] =~ /\d+/){
			$ncores = $ARGV[$i+1];
		}
		else{
			die "Invalid entry '$ARGV[$i+1]' for number of processing cores (-nc). It must be an integer. Please, reset and try again.\n";
		}
	}
	elsif($ARGV[$i] eq "-redo"){
		if($ARGV[$i+1] =~ /^F(ALSE)?/){
			$redo = "";
		}
		elsif($ARGV[$i+1] !~ /^T(RUE)?/){
			die "Invalid entry '$ARGV[$i+1]' for redoing mode (-redo). It must be a boolean (TRUE/FALSE, or shortly T/F). Please, reset and try again.\n";
		}
	}
	else{
		die "Unrecognized flag '$ARGV[$i]': Check for possible typos or missing values.\n\nAll allowed flags and corresponding values:\n\t-s <sequence_alignment_file_name> (fasta or phylip)\n\t-m <substitution_model> (iqtree format)\n\t-t <trees_file> (newick format - one per line)\n\t-n <number_of_replicates> (default = 1000)\n\t-t <tests_demanded> (default = BP/KH/SH)\n\nFor more info, please check the manual.\n";
	}
}
if(!defined($seqs)||!defined($submodel)||!defined($treez)){
	print "Mandatory argument(s) missing:\n";
	if(!defined($seqs)){
		print "\t(-s) sequence alignment file name\n";
	}
	if(!defined($submodel)){
		print "\t(-m) substitution model\n";
	}
	if(!defined($treez)){
		print "\t(-z) trees file name\n";
	}
	die;
}
if(!defined($nreps)){
	print "\tNumber of replicates (-n) not specified. Setting to default '1000'\n";
}
if(!@testconf){
	print "\tTopology tests (-t) not specified. Setting to default 'BP/KH/SH'\n"
}



#armazena os testes requisitados e os procedimentos pelos quais devem ser realizados (aproximação normal, bootstrap paramétrico ou não-paramétrico e níveis de otimização dos parâmetros)
my %ts;
my %ps;
foreach(@testconf){
	my @temp = split(/:/);
	if(exists($temp[1])){
		my @p = sort(split(/,/,$temp[1]));
		foreach(@p){
			$ts{$temp[0]}{$_} = "";
			$ps{$_}{$temp[0]} = "";
			if(($_ == -3 && $temp[0] ne "KH") || ($_>-3 && $_<0 && $temp[0] ne "SOWH") || ($_>=0 && $_ <=2 && $temp[0] eq "SOWH") || ($_>2 && ($temp[0] ne "AU" || $temp[0] ne "BP"))){
				die "Procedure parameter '$_' is not compatible with $temp[0] test.\nPlease, consult manual for a valid procedure, reset required tests argument (-t) and try again.\n";
			}
		}
	}
	elsif($temp[0] eq "SOWH"){
		$ts{SOWH}{-1} = "";
		$ps{-1}{SOWH} = "";
	}
	else{
		$ts{$temp[0]}{0} = "";
		$ps{0}{$temp[0]} = "";
	}
}



#recebe modelo de substituição a ser utilizado para calcular todas as verossimilhanças
if($submodel=~/\{.+\}/){
	die "This version of phyTest cannot handle pre-specified substitution parameter values (any between '{}') in '$submodel'.\nPlease, reset the model argument (-m) and try again.\n";
}
my @model = split(/\+/,$submodel);
my $matrix = shift @model;
my $inv = "";
my $gama = "";
my $free = "";
foreach(@model){
	if(/I/){
		$inv = "+$_";
	}
	elsif(/G\d*/){
		$gama = "+$_";
	}
	elsif(/R\d*/){
		$free = "+$_";
	}
}
if($gama && $free){
	die "Cannot use '+G' and '+R' simultaniously as passed in '$submodel'.\nPlease, reset the model argument (-m) and try again.\n";
}
elsif($free && exists($ts{SOWH})){
	die "Due to limitations in Seq-Gen program, freely distributed rate categories (+R) - as passed in '$submodel' - are not a valid parameter for SOWH's sequence simulations.\nPlease, reset the model (-m) or required tests (-t) and try again.\n";
}



#recebe nome do alinhamento, sem extensão
my $data_name;
if($seqs=~/(.+)\..+/){ # CORRIGIR REGEX PARA TORNAR \..+ OPCIONAL SEM TORNAR O GRUPO DE CAPTURA GREEDY
	$data_name = $1;
}
my $nsites;
my %data;
print "\nReading sequence alignment from '$seqs'...\n";
#se for necessário reamostrar sitios do alinhamento (dado em $seqs) para bootstrap não-paramétrico, identifica seu formato (phylip ou fasta) e o armazena em %data, além do seu número de sítios em $nsites
if(exists($ps{1})||exists($ps{2})||exists($ps{3})||exists($ps{4})){
	open IN, "<$seqs" or die "Failed to open sequence alignment '$seqs': $!";
	my $header;
	my $line1 = <IN>;
	until($line1 !~ /^(\h+)?\n/){
		$line1 = <IN>;
	}
	if($line1 =~ /\d+\h+\d+/){ #se for phylip...
		print "\tPhylip format detected.\n";
		while(<IN>){
			if(/(\w+)\h+([\w\?\*-]+)/){
				$header = $1;
				my @seq = split("",$2);
				$data{$header} = [@seq];
			}
		}
	}
	elsif($line1 =~ /^>\w+/){ #se for fasta...
		print "\tFasta format detected\n";
		seek(IN,-length($line1),1);
		while(<IN>){
			if(/^>(\w+)/){
				$header = $1;
			}
			elsif(/([\w\?\*-]+)/){
				my @seq = split("",$1);
				$data{$header} = [@seq];
			}
		}
	}
	else{ #se não for nenhum dos dois...
		die "Failed to recognize sequence alignment.\nPlease, make sure it is in supported format (either sequential fasta or sequential phylip).\n";
	}
	$nsites = scalar(@{$data{$header}});
	print "\t$nsites sites long\n";
	($header,$line1) = undef;
	close IN;
}
# caso contrário, se for realizar bootstrap paramétrico ou RELL, obtém apenas o número de sítios pelo tamanho da primeira sequencia
elsif(exists($ps{0})||exists($ps{-1})||exists($ps{-2})){
	open IN, "<$seqs";
	while(<IN>){
		if(/^(\w+\h+)?([\w\?\*-]+)/){
			if(defined($1)){
				print "\tPhylip format detected\n";
			}
			else{
				print "\tFasta format detected\n";
			}
			$nsites = length($2);
			last;
		}
	}
	print "\t$nsites sites long\n";
	close IN;
}


#recebe as árvores a serem testadas (em formato newick e listadas em um mesmo doc de texto), associa cada uma a uma ID e as desenraiza
my @trees;
open IN, "<$treez" or die "Failed to open '$treez': $!";
print "\nReading phylogenetic trees from '$treez'...\n";
my $id;
while(<IN>){
	if(/[\(\)\w\.,:]+;/){
		$id++;
		print "\tt${id} = $_\n";
		push(@trees,unroot($_));
	}
}
close IN;
# garante que as árvores recebidas sejam comparáveis, topologicamente únicas e que o número de topologias únicas seja suficiente para os testes requisitados; armazena as árvores em $opt{d}, gerando um arquivo separado para cada uma
my %opt;
if(scalar(@trees)==0){
	die "\nNot enough trees to proceed: no phylogenetic tree detected.\n";
}
else{
	print "\t".scalar(@trees)." tree(s) detected and successfully unrooted. ";
	if(scalar(@trees)>1){
		for($a=0; $a<scalar(@trees); $a++){
			if(defined($trees[$a])){
				$id = "t".($a+1);
				open OUT, ">${id}.tre";
				print OUT $trees[$a];
				close OUT;
				$opt{d}{$id}{tree} = $trees[$a];
				for($b=$a+1; $b<scalar(@trees); $b++){
					if(defined($trees[$b]) && phyDist($trees[$a],$trees[$b])==0){ # ao buscar diferenças topológicas entre as árvores, a função phyDist também garante que sejam comparáveis (mesmo set de táxons)
						$trees[$b] = undef;
					}
				}
			}
		}
		if(scalar(@trees)==scalar(keys(%{$opt{d}}))){
			print "All topologically unique."
		}
		else{
			print "Some have the same topology.\n\tOnly the following ".scalar(keys(%{$opt{d}}))." can be considered: ".join(", ",sort(keys(%{$opt{d}})))."\n";
		}
	}
	if(scalar(@trees)==1 || (exists($opt{d}) && scalar(keys(%{$opt{d}}))==1)){
	 	if(scalar(keys(%ps))>2 || (!exists($ps{3}) && !exists($ps{4}))){
			die "\nNot enough trees to proceed: using the tests and/or procedures passed via (-t), at least 2 reasonable and topologically unique trees must be provided to achieve minimal precision in confidence calculation.\n";
		}
		elsif(scalar(@trees)==1){
			open OUT, ">t1.tre";
			print OUT $trees[0];
			close OUT;
			$opt{d}{t1}{tree} = $trees[0];
		}
	}	
}
undef @trees;

#foreach(2..5,7..10,12..15){
#	delete $opt{d}{"t$_"};
#}

#dado o alinhamento original, utiliza o iqtree para otimizar os comprimentos de ramo, além dos parâmetros de substituição para cada árvore e os armazena junto a suas verossimilhanças em %opt
my $dmlt;
foreach my $t(sort(keys(%{$opt{d}}))){ # para cada árvore $_...
	print "\n\t$t ";
	#via iqtree, otimiza seus comprimentos de ramo e parametros de substituição para o alinhamento dado
	-f "${t}_${data_name}.sitelh" && !$redo ? warn "\tReusing previously computed likelihood of '$t' for sequence alignment '$seqs'.\n" : system "iqtree-omp -s $seqs -te ${t}.tre -m $submodel -pre ${t}_${data_name} -nt 2 -wsl -quiet $redo";

	#rearmazena árvore com comprimentos de ramo otimizados
	open IN, "<${t}_${data_name}.treefile" or die "Failed to open '${t}_${data_name}.treefile': $!";
	$opt{d}{$t}{tree} = <IN>;
	$opt{d}{$t}{tree} = unroot($opt{d}{$t}{tree});
	close IN;

	if(exists($ps{-1})||exists($ps{-2})||exists($ps{1})||exists($ps{3})){ # se for realizar bootstrap paramétrico de qualquer tipo ou não-paramétrico com otimização parcial...
		#armazena parâmetros de susbtituição estimados
		my @rates;
		my @freqs;
		$opt{d}{$t}{pinv} = "";
		$opt{d}{$t}{alpha} = "";
		$opt{d}{$t}{spr} = "";
		open IN, "<${t}_${data_name}.iqtree" or die "Failed to open '${t}_${data_name}.iqtree': $!";
		while(<IN>){
			chomp($_);
			if(/[ACTG]-[ACTG]: (\d\.\d+)/){
				push(@rates,$1);
			}
			elsif(/pi\([ACTG]\) = (\d\.\d+)/){
				push(@freqs,$1);
			}
			elsif(/Proportion of invariable sites: (\d\.\d+)/){
				$opt{d}{$t}{pinv} = $1;
			}
			elsif(/Gamma shape alpha: (\d\.\d+)/){
				$opt{d}{$t}{alpha} = $1;
				last;
			}
			elsif(/Site proportion and rates:\h+(.+)/){
				my $temp = $1;
				$temp =~ s/[\(\)]//;
				$temp =~ s/,/ /;
				my @spr = split(" ",$temp);
				$opt{d}{$t}{spr} = join(",",multiround([@spr],3));
				last;
			}
		}
		close IN;
		if(exists($ps{-1})||exists($ps{1})||exists($ps{3})){ # se realizar otimização parcial...
			#armazena rates no formato utilizado no iqtree
			$opt{d}{$t}{iqrate} = join(",",uniq(grep($_!=1,@rates)));
		}
		if(exists($ps{-2})||exists($ps{-1})){ # se realizar bootstrap paramétrico...
			#armazena rates no formato para restrição do modelo GTR no seq-gen
			$opt{d}{$t}{sgrate} = join(",",@rates);
		}
		$opt{d}{$t}{freq} = join(",",@freqs);
		
	}

	#armazena valor de log-verossimilhança total da árvore e suas log-verossimilhanças por sítio
	open IN, "<${t}_${data_name}.sitelh" or die "Failed to open '${t}_${data_name}.sitelh': $!";
	while(my $line = <IN>){
		if($line =~ /^Site_Lh\h+(.+)/){
			$opt{d}{$t}{slik} = [split(/\h/,$1)];
			$opt{d}{$t}{lik} = round(sum(@{$opt{d}{$t}{slik}}),4);
			last;
		}
	}
	print "lnL = $opt{d}{$t}{lik}";
	# identifica árvore de maior log-verossimilhança para os dados ($dmlt)
	if(!defined($dmlt) || $opt{d}{$dmlt}{lik} < $opt{d}{$t}{lik}){
		$dmlt = $t;
	}
	close IN;

}
print "\n\tHighest-scoring tree: $dmlt (when performing KH, SH or SOWH test, the remaining trees are compared against this)\n";
if(exists($ts{KH})||exists($ts{SH})||exists($ts{SOWH})){ # se for realizado qualquer teste baseado em delta de verossimilhança...
	# armazena diferenças de log-verossimilhaca entre a $dmlt e demais árvores
	foreach my $t(keys(%{$opt{d}})){
		$opt{d}{$t}{delta} = round($opt{d}{$dmlt}{lik}-$opt{d}{$t}{lik},4);
	}
}
my $Coptfreq;
if(exists($ps{-2})||exists($ps{2})|exists($ps{4})){ # se for realizar qualquer tipo de bootstrap com otimização completa...
	if(scalar(uniq(split(",",$opt{d}{t1}{freq})))==1){ # e se as frequencias de bases forem iguais (sinalizando que são restritas no modelo de substituição especificado)...
		# especifica flag para que as frequencias de base sejam fixas (e iguais) durante a reotimização de parametros no iqtree
		$Coptfreq = "+FQ";
	}
	else{ # caso contrário...
		# flag determina que frequencias sejam otimizadas por ML
		$Coptfreq = "+FO";
	}
}




my %pvs;
if(exists($ps{0})||exists($ps{1})||exists($ps{2})||exists($ps{3})||exists($ps{4})){ # se o(s) teste(s) envolver(em) geração de réplicas não-paramétricas...

	my @procs = grep($_>=0,keys(%ps));
	# define uma lista de fatores de escala, para caso venha a ser realizado o teste AU, ou apenas um fator redundante (1), caso não
	my @sf; # fatores de escala
	my @rsf; # raiz quadrada dos fatores de escala
	my @inv_rsf; # inverso da raiz dos fatores de escala
	if(exists($ts{AU})){ 
		@sf = (0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4);
		foreach(@sf){
			my $rf = sqrt($_);
			push(@rsf,$rf);
			push(@inv_rsf,(1/$rf));
		}
	}
	else{
		@sf = (1);
	}
	for(my $k=0; $k<scalar(@sf); $k++){ # para cada fator de escala...

		print "\nGenerating $nreps non-parametric replicates with ".int($nsites*$sf[$k])." sites each and executing any required optimization...\n";
		for(my $r=1; $r<=$nreps; $r++){ # dada cada réplica $r...

			#sorteia, com reposição, as posições a serem amostradas para compor a réplica $r, totalizando N' posições = N posições do alinhamento original ($nsites) X o fator de escala ($sf[$k])
			my @posits;
			for(my $i=0; $i<($nsites*$sf[$k]); $i++){
				$posits[$i] = int(rand($nsites));
			}
			#gera alinhamento-réplica, se necessário para algum dos procedimentos
			if(exists($ps{1})||exists($ps{2})||exists($ps{3})||exists($ps{4})){
				open OUT, ">${data_name}_NPboot${r}.fas";
				foreach my $header(sort(keys(%data))){
					print OUT ">".$header."\n".join("",@{$data{$header}}[@posits])."\n";
				}
				close OUT;
			}

			foreach(@procs){ # para cada procedimento $_...

				if($sf[$k] == 1 || exists($ps{$_}{AU})){ # se o fator de escala atual for cabível,
					my $rmlt = undef;

					if($_ < 3){ # e se o procedimento não envolver busca pela árvore ML,					
						foreach my $t(keys(%{$opt{d}})){ # para cada árvore $t...

							#RELL
							if($_ == 0){
								print "\tResampling estimated log-likelihoods (RELL) for each tree: replicate $r (".$nsites*$sf[$k]." sites long)\n" if $t eq "t1" && $verb;
								#reamostra os valores de verossimilhanças por sítio para cada árvore e armazena os totais resultantes
								$opt{r}{$t}{lik} = round(sum(@{$opt{d}{$t}{slik}}[@posits]),4);

							}
							#Bootstrap não-paramétrico
							else{

								# utiliza o alinhamento-réplica para otimizar os parâmetros de cada árvore
								if($_ == 1){ # com otimização parcial (apenas comprimentos de ramos)...
									print "\tOptimizing trees using replicate $r (".$nsites*$sf[$k]." sites long)\n" if $t eq "t1" && $verb;
									-f "${t}_${data_name}_NPboot${r}_Popt.sitelh" && !$redo ? warn "\tReusing previously computed likelihood of '$t' for a replicate '$r'.\n" : system "iqtree-omp -s ${data_name}_NPboot${r}.fas -te ${t}.tre -m $matrix'{'$opt{d}{$t}{iqrate}'}'$inv'{'$opt{d}{$t}{pinv}'}'$gama'{'$opt{d}{$t}{alpha}'}'$free'{'$opt{d}{$t}{spr}'}'+F'{'$opt{d}{$t}{freq}'}' -pre ${t}_${data_name}_NPboot${r}_Popt -nt 2 -wsl -quiet $redo";
									open IN, "<${t}_${data_name}_NPboot${r}_Popt.sitelh" or die "Failed to open '${t}_${data_name}_NPboot${r}_Popt.sitelh': $!";
								}
								else{ # ou completa (ramos + parametros de susbtituição)
									print "\tOptimizing trees and substitution parameters using replicate $r (".$nsites*$sf[$k]." sites long)\n" if $t eq "t1" && $verb;
									-f "${t}_${data_name}_NPboot${r}_Copt.sitelh" && !$redo ? warn "\tReusing previously computed likelihood of '$t' for a replicate '$r'.\n" : system "iqtree-omp -s ${data_name}_NPboot${r}.fas -te ${t}.tre -m ${matrix}${inv}${gama}${free}${Coptfreq} -pre ${t}_${data_name}_NPboot${r}_Copt -nt 2 -wsl -quiet $redo";
									open IN, "<${t}_${data_name}_NPboot${r}_Copt.sitelh" or die "Failed to open '${t}_${data_name}_NPboot${r}_Copt.sitelh': $!";
								}
								# e armazena sua verossimilhança
								while(my $line = <IN>){
									if($line =~ /^Site_Lh\h+(.+)/){
										$opt{r}{$t}{lik} = round(sum(split(/\h/,$1)),4);
										last;
									}
								}
								close IN;

							}
							# identifica a árvore de maior verossimilhança para a réplica
							if(!defined($rmlt) || $opt{r}{$rmlt}{lik} < $opt{r}{$t}{lik}){
								$rmlt = $t;
							}

						}
					}
					# PARA O BOOTSTRAP COM BUSCA TOPOLÓGICA PODE SER INTERESSANTE IMPLEMENTAR A BUSCA TAMBÉM COM OS DADOS ORIGINAIS E UTILIZAR OS PARAMETROS DE SUBSTITUIÇÃO ESTIMADOS PARA A BUSCA COM OTIMIZAÇÃO PARCIAL A SEGUIR
					else{ # se o procedimento envolver busca pela árvore ML,
						# utiliza o alinhamento-réplica para inferí-la
						if($_ == 3){ # com otimização parcial (apenas topologia + comprimentos de ramos)
							print "\tSearching ML tree for replicate $r (".$nsites*$sf[$k]." sites long)\n" && $verb;
							-f "MLt_${data_name}_NPboot${r}_TrSch_Popt.sitelh" && !$redo ? warn "\tReusing previously computed maximum likelihood for a replicate '$r'.\n" : system "iqtree-omp -s ${data_name}_NPboot${r}.fas -m $matrix'{'$opt{d}{$dmlt}{iqrate}'}'$inv'{'$opt{d}{$dmlt}{pinv}'}'$gama'{'$opt{d}{$dmlt}{alpha}'}'$free'{'$opt{d}{$dmlt}{spr}'}'+F'{'$opt{d}{$_}{freq}'}' -pre MLt_${data_name}_NPboot${r}_TrSch_Popt -nt 2 -wsl -quiet $redo";
							open IN, "<MLt_${data_name}_NPboot${r}_TrSch_Popt.treefile" or die "Failed to open 'MLT_${data_name}_NPboot${r}_TrSch_Popt.sitelh': $!";
							$rmlt = <IN>;
							close IN;
						}
						elsif($_ == 4){ # ou completa (topologia + ramos + parametros de substituição)
							print "\tSearching ML tree and substitution parameters for replicate $r (".$nsites*$sf[$k]." sites long)\n" && $verb;
							-f "MLt_${data_name}_NPboot${r}_TrSch_Copt.sitelh" && !$redo ? warn "\tReusing previously computed maximum likelihood for a replicate '$r'.\n" : system "iqtree-omp -s ${data_name}_NPboot${r}.fas -m ${matrix}${inv}${gama}${free}${Coptfreq} -pre MLt_${data_name}_NPboot${r}_TrSch_Copt -nt 2 -wsl -quiet $redo";
							open IN, "<MLt_${data_name}_NPboot${r}_TrSch_Copt.treefile" or die "Failed to open 'MLT_${data_name}_NPboot${r}_TrSch_Copt.sitelh': $!";
							$rmlt = <IN>;
							close IN;
						}
					}

					my $elw_ratiosum = 0 if $sf[$k] == 1 && exists($ps{$_}{ELW});
					foreach my $t(keys(%{$opt{d}})){

						if($sf[$k] == 1){
							if(exists($ps{$_}{KH})){ # se foi requisitado KH sob o procedimento
								push(@{$pvs{$t}{KH}{$_}{rdeltas}},$opt{r}{$dmlt}{lik}-$opt{r}{$t}{lik}); # armazena delta entre a árvore que foi a de ML para o alinhamento original e $t (para posterior centralização)
							}
							if(exists($ps{$_}{SH})){ # se foi requisitado SH
								push(@{$pvs{$t}{SH}{$_}{rliks}},$opt{r}{$t}{lik}); # armazena apenas lnL de $t (para posterior centralização)
							}
							if(exists($ps{$_}{ELW})){ # se foi requisitado ELW
								$pvs{$t}{ELW}{$_}{rratio} = exp($opt{r}{$t}{lik}-$opt{r}{$rmlt}{lik}); # armazena razão entre a lnL de $t e da árvore de ML para $r (para computação mais precisa da razão, é obtida de e^delta, sendo delta = diferença entre lnLs)
								$elw_ratiosum += $pvs{$t}{ELW}{$_}{rratio}; # e a acumula à soma das razões para as diferentes árvores
							}
						}
						if(exists($ps{$_}{BP})||exists($ps{$_}{AU})){ # se foi requisitado BP ou AU sob o procedimento,
							$pvs{$t}{BP}{$_}{$k} = 0 if !exists($pvs{$t}{BP}{$_}{$k});
							if($_ < 3 && $t eq $rmlt){ # sendo ele RELL, otimização parcial ou total e se $t tiver obtido maior lnL para a réplica atual,
								$pvs{$t}{BP}{$_}{$k}++; # adiciona à contagem de vezes em que $t foi a árvore de ML
							}
							elsif($_ >= 3 && phyDist($opt{d}{$t}{tree},$rmlt) == 0){ # ou sendo ele busca topológica + otimização parcial ou total e se $t tiver a topologia de ML para a réplica atual,
								$pvs{$t}{BP}{$_}{$k}++; # adiciona à contagem de vezes em que $t foi a árvore de ML
							}
						}

					}
					# terminadas todas as árvores
					if($sf[$k] == 1 && exists($ps{$_}{ELW})){ # se foi requisitado ELW sob o procedimento
						foreach my $t(keys(%{$opt{d}})){ # novamente, dada cada árvore $t...
							$pvs{$t}{ELW}{$_}{weightsum} = 0 if !exists($pvs{$t}{ELW}{$_}{weightsum});
							$pvs{$t}{ELW}{$_}{weightsum} += round(round($pvs{$t}{ELW}{$_}{rratio},4)/round($elw_ratiosum,4),4); # armazena peso da razão de sua verossimilhança (ponderado pela soma das razões, calculada acima)
							delete $pvs{$t}{ELW}{$_}{rratio};
						}
						$elw_ratiosum = undef;
					}

				}
			}

		}
		# terminadas todas as réplicas do fator $k, calcula pvalues (de BP, KH, SH e ELW apenas)
		foreach(@procs){ # novamente, para cada procedimento $_...

			foreach my $t(keys(%{$opt{d}})){ # e para cada árvore $t...

				if($sf[$k] == 1){
					if(exists($ps{$_}{BP})){ # se foi requisitado BP

						print "\tComputing confidence of each tree via BP:$_\n" if $t eq "t1" && $verb;
						$pvs{$t}{BP}{$_}{pval} = $pvs{$t}{BP}{$_}{$k}/$nreps;

					}
					if(exists($ps{$_}{KH})){ # se foi requisitado KH

						print "\tComputing confidence of each tree via KH:$_\n" if $t eq "t1" && $verb;
						$pvs{$t}{KH}{$_}{pval} = 0 if !exists($pvs{$t}{KH}{$_}{pval});
						my $rdelta_mean = round((sum(@{$pvs{$t}{KH}{$_}{rdeltas}})/$nreps),4);
						for(my $r=1; $r<=$nreps; $r++){ # para cada réplica $r...
							$pvs{$t}{KH}{$_}{rdeltas}[$r] -= $rdelta_mean; # executa o procedimento de centralização de deltas							
							if($opt{d}{$t}{delta}<=$pvs{$t}{KH}{$_}{rdeltas}[$r]){ # se o delta entre a árvore de ML original e $t superar o delta original...
								$pvs{$t}{KH}{$_}{pval}++; # adiciona à contagem de vezes em que o delta original foi superado
							}
						}
						$pvs{$t}{KH}{$_}{pval} = round($pvs{$t}{KH}{$_}{pval}/$nreps,4);
						delete $pvs{$t}{KH}{$_}{rdeltas};

					}
					if(exists($ps{$_}{SH})){ # se for requisitado SH

						$pvs{$t}{SH}{$_}{pval} = 0 if !exists($pvs{$t}{SH}{$_}{pval});
						$pvs{$t}{SH}{$_}{likmean} = round((sum(@{$pvs{$t}{SH}{$_}{rliks}})/$nreps),4); # armazena a média das lnLs entre réplicas temporariamente

					}
					if(exists($ps{$_}{ELW})){ # se foi requisitado ELW

						print "\tComputing confidence of each tree via ELW:$_\n" if $t eq "t1" && $verb;
						$pvs{$t}{ELW}{$_}{pval} = round($pvs{$t}{ELW}{$_}{weightsum}/$nreps,4);
						delete $pvs{$t}{ELW}{$_}{weightsum};

					}
				}
				if(exists($ps{$_}{AU})){ # se foi requisitado AU sob o procedimento

					my $bp = $pvs{$t}{BP}{$_}{$k}/$nreps;
					my ($cq,$w);
					$pvs{$t}{AU}{$_}{kn} = 0 if !exists($pvs{$t}{AU}{$_}{kn});
					if($bp == 0 || $bp == 1){
						$cq = $w = 0;
					}else{
						$cq = inv_cdf(1-$bp); # quantil do complemento do BP
						$w = ((pdf($cq)**2)*$nreps)/($bp*(1-$bp)); #least-squares weight						
						$pvs{$t}{AU}{$_}{kn}++;
					}
					my $wa = $w*$rsf[$k];
					my $wb = $w*$inv_rsf[$k];
					exists($pvs{$t}{AU}{$_}{wab}) ? $pvs{$t}{AU}{$_}{wab} += $wa*$inv_rsf[$k] : $pvs{$t}{AU}{$_}{wab} = $wa*$inv_rsf[$k];
					exists($pvs{$t}{AU}{$_}{wbc}) ? $pvs{$t}{AU}{$_}{wbc} += $wb*$cq : $pvs{$t}{AU}{$_}{wbc} = $wb*$cq;
					exists($pvs{$t}{AU}{$_}{wac}) ? $pvs{$t}{AU}{$_}{wac} += $wa*$cq : $pvs{$t}{AU}{$_}{wac} = $wa*$cq;
					exists($pvs{$t}{AU}{$_}{waa}) ? $pvs{$t}{AU}{$_}{waa} += $wa*$rsf[$k] : $pvs{$t}{AU}{$_}{waa} = $wa*$rsf[$k];
					exists($pvs{$t}{AU}{$_}{wbb}) ? $pvs{$t}{AU}{$_}{wbb} += $wb*$inv_rsf[$k] : $pvs{$t}{AU}{$_}{wbb} = $wb*$inv_rsf[$k];
					#delete $pvs{$t}{BP}{$_}{$k};

				}

			}
			# terminadas todas as árvores
			if($sf[$k] == 1 && exists($ps{$_}{SH})){ # se foi requisitado SH
				
				for(my $r=1; $r<=$nreps; $r++){ # para cada réplica $r...
					my $rmlt;
					foreach my $t(keys(%{$opt{d}})){ # e para cada árvore $t...
						$pvs{$t}{SH}{$_}{rliks}[$r] -= $pvs{$t}{SH}{$_}{likmean}; # executa o procedimento de centralização de lnLs
						# seleciona a árvore que obtém maior lnL centralizado para a réplica
						if(!defined($rmlt) || $pvs{$rmlt}{SH}{$_}{rliks}[$r] < $pvs{$t}{SH}{$_}{rliks}[$r]){
							$rmlt = $t;
						}
					}
					foreach my $t(keys(%{$opt{d}})){ # para cada árvore, novamente...
						if($opt{d}{$t}{delta}<=($pvs{$rmlt}{SH}{$_}{rliks}[$r]-$pvs{$t}{SH}{$_}{rliks}[$r])){ # se o delta entre e a árvore de ML para $r e $t superar o delta original...
							$pvs{$t}{SH}{$_}{pval}++; # adiciona à contagem de vezes em que o delta original foi superado
						}
					}
				}
				print "\tComputing confidence of each tree via SH:$_\n" if $verb;
				foreach my $t(keys(%{$opt{d}})){
					$pvs{$t}{SH}{$_}{pval} = round($pvs{$t}{SH}{$_}{pval}/$nreps,4);
					delete $pvs{$t}{SH}{$_}{likmean};
					delete $pvs{$t}{SH}{$_}{rliks};
				}

			}

		}

	}
	#terminados todos os fatores de escala, calcula pvalue do AU
	if(exists($ts{AU})){

		foreach(@procs){
			print "\tComputing confidence of each tree via AU:$_\n" if $verb;
			foreach my $t (sort(keys(%{$opt{d}}))){
				if($pvs{$t}{AU}{$_}{kn} >= 2){
					my $denom = ($pvs{$t}{AU}{$_}{wab}**2)-($pvs{$t}{AU}{$_}{waa}*$pvs{$t}{AU}{$_}{wbb});
					#$pvs{$t}{AU}{$_}{se} = 0;
					if($denom == 0){
						delete $pvs{$t}{AU}{$_};
						$pvs{$t}{AU}{$_}{pval} = 0;
					}else{
						my $D = (($pvs{$t}{AU}{$_}{wab}*$pvs{$t}{AU}{$_}{wbc})-($pvs{$t}{AU}{$_}{wbb}*$pvs{$t}{AU}{$_}{wac}))/$denom;
						my $C = (($pvs{$t}{AU}{$_}{wab}*$pvs{$t}{AU}{$_}{wac})-($pvs{$t}{AU}{$_}{waa}*$pvs{$t}{AU}{$_}{wbc}))/$denom;
						my $std_err = -($pvs{$t}{AU}{$_}{waa}+$pvs{$t}{AU}{$_}{wbb}+(2*$pvs{$t}{AU}{$_}{wab}))/$denom;
						delete $pvs{$t}{AU}{$_};
						$pvs{$t}{AU}{$_}{pval} = 1-cdf($D-$C);
						#$pvs{$t}{AU}{$_}{se} = pdf($D-$C)*sqrt($std_err);
					}
					#IMPLEMENTAR cálculo e teste de significância da sum of square differences (RSS) - se for muito alta, AU não deve ser utilizado (Shimodaira, 2002)
					# cálculo do RSS segue: sum_k( $w[$k] * ( ($D*$rsf[$k] + $C*$inv_rsf[$k]) - $cq[$k] )**2 )
				}else{
					$pvs{$t}{AU}{$_}{pval} = 0.77;
=cut
					my $cq_sum = 0;
					for(my $k=0;$k<scalar(@sf);$k++){
						my $bp = $pvs{$t}{BP}{$_}{$k}/$nreps;
						$cq_sum += ($bp == 0 || $bp == 1) ? 0 : inv_cdf(1-$bp);
					}
					$pvs{$t}{AU}{$_}{pval} = ($cq_sum >= 0) ? 0 : 1;
=cut
				}
				delete $pvs{$t}{BP}{$_} if !exists($ps{$_}{BP});
				delete $pvs{$t}{BP} if !exists($ts{BP});
			}
		}

	}

}
if(exists($ps{-1})||exists($ps{-2})){ # se envolver(em) geração de réplicas paramétricas...

	my @procs = grep(/-(1|2)/,keys(%ps));
	# converte número de categorias gama para o formato do seq-gen
	my $gcat;
	if(defined($gama)){
		if($gama =~ /\+G(\d+)/){
			$gcat = "-g $1";
		}
		else{
			$gcat = "-g 4";
		}
	}

	foreach my $truet (keys(%{$opt{d}})){ # para cada árvore $truet

		if($truet != $dmlt){ # que não for a de ML para o alinhamento original,

			# converte alpha e proporção de sítios invariáveis
			my $alph;
			if(defined($gama)){
				$alph = "-a $opt{d}{$truet}{alpha}";
			}
			my $pinv;
			if(defined($inv)){
				$pinv = "-i $opt{d}{$truet}{pinv}";
			}

			# e simula N alinhamentos utilizando a árvore $truet, sendo N = número de réplicas dado em $nreps
			print "\nUsing $truet to simulate $nreps parametric replicates with ".$nsites." sites each and reoptimizing the other trees...\n";
			system "seq-gen -m GTR -l $nsites -n $nreps $gcat $alph $pinv -f $opt{d}{$truet}{freq} -r $opt{d}{$truet}{sgrate} -or < ${truet}.tre > ${data_name}-${truet}_Pboots.data";	
			open IN, "<${data_name}-${truet}_Pboots.data" or die "Failed to open '${data_name}-${truet}_Pboots.data': $!";
			local $/ = "\n ";
			my $r;	
			while(<IN> =~ /\d+\h\d+\n(\w+\h[ACGT]+\n)+/){ # para cada alinhamento gerado...
				$r++;

				# cria arquivo separado o contendo
				open OUT, ">${data_name}-${truet}_Pboot${r}.phy";
				print OUT;
				close OUT;
				foreach(@procs){ # para cada procedimento $_...

					my $rmlt;
					foreach my $t(keys(%{$opt{d}})){ # e, para cada árvore, utiliza o alinhamento simulado para otimização

						if($_ == -1){	# parcial (apenas dos comprimentos de ramos, fixando os mesmos parâmetros de substituição utilizados para a simulação dos dados)
							-f "${t}_${data_name}-${truet}_Pboot${r}_Popt.sitelh.sitelh" && !$redo ? warn "\tReusing previously computed likelihood of '$t' for a replicate '$r'.\n" : system "iqtree-omp -s ${data_name}-${truet}_Pboot${r}.phy -te ${t}.tre -m $matrix'{'$opt{d}{$truet}{iqrate}'}'$inv'{'$opt{d}{$truet}{pinv}'}'$gama'{'$opt{d}{$truet}{alpha}'}'+F'{'$opt{d}{$truet}{freq}'}' -pre ${t}_${data_name}-${truet}_Pboot${r}_Popt -nt 2 -wsl -quiet $redo";
							open IN2, "<${t}_${data_name}-${truet}_Pboot${r}_Popt.sitelh" or die "Failed to open '${t}_${data_name}-${truet}_Pboot${r}_Popt.sitelh': $!";
						}
						elsif($_ == -2){	# ou completa (ramos + parametros de susbtituição)
							-f "${t}_${data_name}-${truet}_Pboot${r}_Copt.sitelh" && !$redo ? warn "\tReusing previously computed likelihood of '$t' for a replicate '$r'.\n" : system "iqtree-omp -s ${data_name}-${truet}_Pboot${r}.phy -te ${t}.tre -m ${matrix}${inv}${gama}${Coptfreq} -pre ${t}_${data_name}-${truet}_Pboot${r}_Copt -nt 2 -wsl -quiet $redo";
							open IN2, "<${t}_${data_name}-${truet}_Pboot${r}_Copt.sitelh" or die "Failed to open '${t}_${data_name}-${truet}_Pboot${r}_Copt.sitelh': $!";
						}
						#calcula verossimilhança de 
						while(<IN2>){
							if(/^Site_Lh\h+(.+)/){
								$opt{r}{$t}{lik} = round(sum(split(/\h/,$1)),4);
								last;
							}
						}
						close IN2;
						if(!defined($rmlt) || $opt{0-$_}{$rmlt}{lik} < $opt{r}{$t}{lik}){
							$rmlt = $t;
						}

					}

					$pvs{$truet}{SOWH}{$_}{pval} = 0 if !exists($pvs{$truet}{SOWH}{$_}{pval});
					if($opt{d}{$truet}{lik}{delta}<$opt{r}{$rmlt}{lik}-$opt{r}{$truet}{lik}){ # se o delta entre e a árvore de ML para $r e $truet superar o delta original...
						$pvs{$truet}{SOWH}{$_}{pval}++; # adiciona à contagem de vezes em que o delta original foi superado
					}

				}

			}
			close IN;

			#terminadas todas as réplicas paramétricas, calcula o pvalue de $truet
			print "\tComputing confidence of each tree via SOWH:$_\n";
			foreach(@procs){ # novamente, para cada procedimento $_...
				$pvs{$truet}{SOWH}{$_}{pval} = $pvs{$truet}{SOWH}{$_}{pval}/$nreps;
			}

		}

	}

}
delete $opt{r};
if(exists($ps{-3})){ # se envolver(em) presunção de normalidade de delta...

	print "\tComputing confidence of each tree via KH:-3\n";
	foreach my $t(keys(%{$opt{d}})){ # para cada árvore $t...

		# armazena vetor de deltas por sítio,
		my @sdeltas;
		foreach(my $s=0;$s<$nsites;$s++){
			push(@sdeltas,round($opt{d}{$dmlt}{slik}[$s]-$opt{d}{$t}{slik}[$s]),4);
		}
		# calcula sua média,
		my $sdelta_mean = round(sum(@sdeltas),4)/$nsites;
		# calcula sua variância,
		my $sdelta_var;
		foreach(my $s=0;$s<$nsites;$s++){
			$sdelta_var += ($sdeltas[$s]-$sdelta_mean)**2
		}
		$sdelta_var = round($sdelta_var,4);
		$sdelta_var /= $nsites-1;
		# aproxima o desvio padrão de delta (sob a H0): DP dos deltas por sítio x raiz do número de sítios (Goldman, 2000)
		my $aprx_delta_sd = sqrt($sdelta_var); # CONFERIR SE É sqrt($sdelta_var) OU  sqrt($sdelta_var*$nsites) EM Kihino & Hasegawa (1989)
		# calcula o pvalor do delta original
		$pvs{$t}{KH}{-3}{pval} = round(1-cdf($opt{d}{$t}{delta},0,$aprx_delta_sd),4);

	}

}


#imprime tabela com a log-verossimilhança, o delta em relação à maior verossimilhança e o p-value de cada árvore pelos diferentes testes requisitados
print "\ntree_id\tlnL";
print "\tdelta" if exists($ts{KH}) || exists($ts{SH}) || exists($ts{SOWH});
foreach my $test(sort(keys(%{$pvs{t1}}))){
	foreach(sort(keys(%{$pvs{t1}{$test}}))){
		print "\t$test:$_";
	}
}
foreach my $tree(sort(keys(%pvs))){
	print "\n$tree\t$opt{d}{$tree}{lik}";
	printf("\t%.4f",$opt{d}{$tree}{delta}) if exists($ts{KH}) || exists($ts{SH}) || exists($ts{SOWH});;
	foreach my $test(sort(keys(%{$pvs{$tree}}))){
		foreach(sort(keys(%{$pvs{$tree}{$test}}))){
			printf("\t%.4f",$pvs{$tree}{$test}{$_}{pval});
		}
	}
}
print "\n\n";

# ==================== FIM DO PROGRAMA ==================== #


# ======================== Math::Gauss ======================== #
# A SINTAXE DAS FUNÇÕES FOI MODIFICADA PARA POUPAR ESPAÇO
sub pdf {
 	my $x = shift if @_;
 	my $m = @_ ? shift : 0;
 	my $s = @_ ? shift : 1;
	die("Can't evaluate Math::Gauss:pdf for \$s=$s not strictly positive") if $s <= 0;
 	my $z = ($x-$m)/$s;
 	return exp(-0.5*$z*$z)/(2.506628274631*$s);
}
sub cdf {
 	my $x = shift if @_;
 	my $m = @_ ? shift : 0;
 	my $s = @_ ? shift : 1;
 	die( "Can't evaluate Math::Gauss:cdf for \$s=$s not strictly positive") if $s <= 0;
 	my $z = ($x-$m)/$s;
 	my $t = 1.0/(1.0 + 0.2316419*abs($z));
 	my $y = $t*(0.319381530
	    + $t*(-0.356563782
		  + $t*(1.781477937
			+ $t*(-1.821255978
				+ $t*1.330274429))));
 	return $z > 0 ? 1.0-pdf($z)*$y : pdf($z)*$y;
}
sub inv_cdf { # MODIFICADA PARA RETORNAR QUANTIL TAMBÉM PARA NORMAIS NÃO-PADRÃO, PORÉM DESVIA DO QNORM DO R - RECOMENDADO ARREDONDAMENTO DO QUANTIL ATÉ SEGUNDA CASA
 	my $p = shift if @_;
 	my $m = @_ ? shift : 0;
 	my $s = @_ ? shift : 1;
 	die("Can't evaluate Math::Gauss::inv_cdf for \$p=$p outside ]0,1[") if $p<=0.0 || $p>=1.0;
 	my $t;
 	$t = $p<0.5 ? sqrt(-2.0*log($p)) : sqrt(-2.0*log(1.0-$p));
 	my $y = (2.515517 + $t*(0.802853 + $t*0.010328));
 	$y /= 1.0 + $t*(1.432788 + $t*(0.189269 + $t*0.001308));
 	my $z;
 	$z = $p<0.5 ? $y-$t : $t-$y;
 	return ($z*$s)+$m;
}
# ============================================================ #

# ================ SUBROTINAS DE ARREDONDAMENTO ============== #
sub round { #arredondamento (5 to up) para um float
	my $float = shift @_;
	my $place_lim = @_ ? shift : 3;
	return $float if int($float)-$float == 0;
	my @sides = split(/\./,$float);
	my @temp = split(/e/,$sides[1]);
	my $note = "";
	if(defined($temp[1])){
		$temp[1] > 0 || 0-$temp[1] < $place_lim ? $note = "e".$temp[1] : return "0";
	}
	return $float if $place_lim >= length($temp[0]);
	my @places = split("",$temp[0]);
	for(my $p=$place_lim;$p<scalar(@places);$p++){
		return (($sides[0].".".join("",@places[0..$place_lim-1]))+(1/(10**$place_lim))).$note if $places[$p]>4;
		return $sides[0].".".join("",@places[0..$place_lim-1]).$note if $places[$p]<4;
	}
	return $sides[0].".".join("",@places[0..$place_lim-1]).$note;
}
sub multiround { #arredondamento (5 to up) para múltiplos floats
	my @rounded;
	my $place_lim = defined($_[1]) ? pop : 3;
	foreach(@{$_[0]}){
		do {(push(@rounded,$_));last} if int($_)-$_ == 0;
		my @sides = split(/\./,$_);
		my @temp = split(/e/,$sides[1]);
		my $note = "";
		if(defined($temp[1])){
			$temp[1] >= 0 || $place_lim >= 0-$temp[1] ? $note = "e".$temp[1] : push (@rounded,0);
		}
		do {(push(@rounded,$_));last} if $place_lim >= length($temp[0]);
		my @places = split("",$temp[0]);
		for(my $p=$place_lim;$p<scalar(@places);$p++){
			do {push(@rounded,(($sides[0].".".join("",@places[0..$place_lim-1]))+(1/(10**$place_lim))).$note);last} if $places[$p]>4;
			do {push(@rounded,($sides[0].".".join("",@places[0..$place_lim-1])).$note);last} if $places[$p]<4 || ($p+1)==scalar(@places);
		}
	}
	return @rounded;
}
# ============================================================ #

# ================== SUBROTINAS TOPOLÓGICAS ================== #
sub unroot { ## Desenraiza árvore em formato newick!
  chomp($_[0]);
  my $t = $_[0];
  substr($t,0,1) = "";
  substr($t,-2,1) = "";
  my $sep = "";
  $t =~ /\:/ ? $sep = ":" : $t =~ s/(\w+|\))/$1:/g;
  my @tree = split(":",$t);
  my (@subtrees,@edges,$st);
  my $o = 0;
  for(my $i=0; $i<scalar(@tree)-1; $i++){
    $st .= "$tree[$i]$sep";
    if($tree[$i] =~ /(\(+)/){
      $o += length($1);
    }
    elsif($tree[$i] =~ /\)+/){
      $o--;
    }
    if($o==0){
      $st =~ s/:$//;
      push(@subtrees,$st);
      my @node = split(/,|;/,$tree[$i+1]);
      push(@edges,$node[0]);
      $tree[$i+1] = $node[1] if defined($node[1]);
      $st = undef;
    }
  }
  if(scalar(@subtrees)==2){
    for(my $i=0;$i<=1;$i++){
      if($subtrees[$i] =~ /\(.+\)/){
        $subtrees[$i] = substr($subtrees[$i],1,length($subtrees[$i])-2);
        do {my $newedge = $edges[0]+$edges[1];$subtrees[$i-1] .= "$sep$newedge"} if $sep eq ":";
        return("(".join(",",@subtrees).");");
      }
    }
  }
  elsif(scalar(@subtrees)>2){
    return($_[0]);
  }
  else{
    die("Tree '$_[0]' is malformed; unable to unroot it.\n");
  }
}
sub phyDist{ ## Dá a distância topológica entre duas árvores newick não enraizadas, pela métrica de partição (sensível apenas a topologia) ou BSD (sensível também a comprimentos de ramo)
	my $lensense = scalar(@_)==3 ? pop : 0; ## 0 - Partition Metric (Penny & Hendy, 1985) =default= | 1 - Branch Score Distance (Kuhner & Felsenstein, 1994)
	my %phys;
	my $i;
	foreach my $t(@_){
		$i++;
		chomp($t);
		my $bl = 0;
		$t =~ /\:/ ? $bl = 1 : $t =~ s/(\w+|\))/$1:/g;
		die("Trees '$_[0]' and '$_[1]' are not comparable: BSD metric not valid for lengthless trees (try default option 0 - Penny & Handy's Partition Metric)") if $bl==0 && $lensense==1;
		my @tree = split(":",$t);
		my ($c,$cc) = 0;
		my (@taxa,@term_bls,@clades,@int_bls,@ocs);
		foreach my $slice(@tree){
			if($slice =~ /(\d\.?\d*)?(\))?,?(\(+)?(\w+)?/){
				if(defined($1)){
					!defined($cc) || $cc == $ocs[-1] ? push(@term_bls,$1) : do {$int_bls[$cc] = $1;$cc = $ocs[-1]};
				}
				if(defined($4)){
					push(@taxa,$4);
					do {$cc = $c+length($3)-1;push(@ocs,($c..$cc));@clades[$c..$cc] = [()];$c+=length($3)} if defined($3);
					foreach(@ocs){
						push(@{$clades[$_]},$4);
					}
				}
				elsif(defined($2)){
					pop(@ocs);
					$clades[$cc] = join(",",sort(@{$clades[$cc]}));
					$cc = $ocs[-1] if $bl == 0;
				}
			}		
		}
		shift(@clades);
		shift(@int_bls);
		for(my $j=0;$j<scalar(@taxa);$j++){
			$phys{$i}{$taxa[$j]} = shift @term_bls;
		}
		for(my $j=0;$j<scalar(@clades);$j++){
			$phys{$i}{$clades[$j]} = shift @int_bls;
		}
	}
	my $PD = 0;
	my @t = (1,2);
	for(my $i=0; $i<2; $i++){
		foreach(keys(%{$phys{$t[$i]}})){
			if(exists($phys{$t[$i-1]}{$_})){
				do {$PD += ($phys{$t[$i]}{$_}-$phys{$t[$i-1]}{$_})**2;$phys{$t[$i]}{$_} = 0;$phys{$t[$i-1]}{$_} = 0} if $lensense == 1;
			}
			else{
				die("Trees '$_[0]' and '$_[1]' are not comparable: different sets of taxa.\n") if $_ !~ /,/;
				$PD += $lensense == 1 ? $phys{$t[$i]}{$_}**2 : 1;
			}
		}
	}
	return($PD);
}
# ============================================================ #

exit;
