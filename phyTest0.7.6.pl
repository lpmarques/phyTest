#!/usr/bin/perl

=cut

Title: phyTest

Version: 0.7 (alpha)

	0.7.1: few modulatizations
	0.7.2: additional modularizations and inclusion of a full tree search to determine subst. parameters along with procedure 3
	0.7.3: PT now reveals seed and allows for user provided seed
	0.7.4: more efficient switches to id required tests and procedures
	0.7.5: natural sorting of tree ids
	0.7.6: corrections to phyDist function and to procedures 3 and 4

Author: Lucas Marques

=cut

use strict;
use warnings;
use List::Util "sum";
use List::MoreUtils "uniq";
use Data::Dumper "Dumper";

# ====================== Math::Gauss ====================== #
# A SINTAXE DAS FUNÇÕES FOI MODIFICADA PARA POUPAR ESPAÇO
sub pdf{
 	my $x = shift;
 	my $m = @_ ? shift : 0;
 	my $s = @_ ? shift : 1;
	die("Can't evaluate Math::Gauss:pdf for \$s=$s not strictly positive") if $s <= 0;
 	my $z = ($x-$m)/$s;
 	return exp(-0.5*$z*$z)/(2.506628274631*$s);
}

sub cdf{
 	my $x = shift;
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

sub inv_cdf{ # MODIFICADA PARA RETORNAR QUANTIL TAMBÉM PARA NORMAIS NÃO-PADRÃO, PORÉM DESVIA DO QNORM DO R - RECOMENDADO ARREDONDAMENTO DO QUANTIL ATÉ SEGUNDA CASA
 	my $p = shift;
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

# =================== ROUNDING SUBROUTINES =================== #
sub round{ #arredondamento (5 to up) para um float
	my $float = shift @_;
	my $place_lim = @_ ? shift : 3;
	return $float if int($float)-$float == 0;
	my @sides = split(/\./,$float);
	return $sides[0] if !defined($sides[1]);
	my @temp = split(/e/,$sides[1]);
	my $note = "";
	if(defined($temp[1])){
		if($temp[1] > 0){
			return $float;
		}
		elsif(0-$temp[1] < $place_lim){
			return 0;
		}
		else{
			$note = "e".$temp[1];	
		}
	}
	return $float if $place_lim >= length($temp[0]);
	my @places = split("",$temp[0]);
	for(my $p=$place_lim;$p<scalar(@places);$p++){
		return (($sides[0].".".join("",@places[0..$place_lim-1]))+(1/(10**$place_lim))).$note if $places[$p]>4;
		return $sides[0].".".join("",@places[0..$place_lim-1]).$note if $places[$p]<4;
	}
	return $sides[0].".".join("",@places[0..$place_lim-1]).$note;
}

sub multiround{ #arredondamento (5 to up) para múltiplos floats
	my @rounded;
	my $place_lim = defined($_[1]) ? pop : 3;
	foreach(@{$_[0]}){
		do {(push(@rounded,$_));last} if int($_)-$_ == 0;
		my @sides = split(/\./,$_);
		do {(push(@rounded,$_));last} if !defined($sides[1]);
		my @temp = split(/e/,$sides[1]);
		my $note = "";
		if(defined($temp[1])){
			if($temp[1] > 0){
				push(@rounded,$_);
			}
			elsif(0-$temp[1] < $place_lim){
				push(@rounded,0);
			}
			else{
				$note = "e".$temp[1];	
			}
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

# ================== TOPOLOGIC SUBROUTINES =================== #
sub unroot{ ## Desenraiza árvore em formato newick!
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
	my $lensense = defined($_[2]) ? pop : 0; ## 0 - Partition Metric (Penny & Hendy, 1985 [or Robindon & Foulds, 1981, originally]) =default= | 1 - Branch Score Distance (Kuhner & Felsenstein, 1994)
	my $i;
	my %parts;
	foreach my $t(@_){
		$i++;
		chomp($t);
		my $bl = 0;
		$t =~ /\:/ ? $bl = 1 : $t =~ s/(\w+|\))/$1:/g;
		die("Trees '$_[0]' and '$_[1]' are not comparable: BSD metric not valid for lengthless trees (try default option 0 - Penny & Handy's Partition Metric)") if !$bl && $lensense;
		my @tree = split(":",$t);
		my ($c,$cci) = 0; # counter of clades and current clade index 
		my (@taxa,@term_bls,@clades,@int_bls,@ocis); # taxa, terminal branch lengths, clades, internal branch lengths and open clades index
		foreach my $slice(@tree){
			if($slice =~ /(\d\.?\d*)?(\))?,?(\(+)?(\w+)?/){
				if(defined($1)){
					!defined($cci) || $cci == $ocis[-1] ? push(@term_bls,$1) : do {$int_bls[$cci] = $1;$cci = $ocis[-1]};
				}
				if(defined($4)){
					push(@taxa,$4);
					do {$cci = $c+length($3)-1;push(@ocis,($c..$cci));@clades[$c..$cci] = [()];$c+=length($3)} if defined($3);
					foreach(@ocis){
						push(@{$clades[$_]},$4);
					}
				}
				elsif(defined($2)){
					pop(@ocis);
					$clades[$cci] = [sort(@{$clades[$cci]})]; 
					$cci = $ocis[-1] if !$bl;
				}
			}
		}
		shift(@clades);
		shift(@int_bls);

		my ($j,$k,$l);
		my $ntaxa = scalar(@taxa);
		my $nclades = scalar(@clades);
		for($j=0;$j<$nclades;$j++){
			my $cntaxa = scalar(@{$clades[$j]});
			if($cntaxa>$ntaxa/2){
				my @temp;
				for($k=0;$k<$ntaxa;$k++){
					for($l=0;$l<$cntaxa;$l++){
						last if $taxa[$k] eq $clades[$j][$l];
					}
					push(@temp,$taxa[$k]) if $l==$cntaxa;
				}
				$clades[$j] = join(",",sort(@temp));
			}
			else{
				$clades[$j] = join(",",@{$clades[$j]});
			}
			$parts{"t$i"}{$clades[$j]} = $bl ? $int_bls[$j] : 0; # stores each internal branch partition (id: list of taxon names in smaller subtree)

		}
		for($j=0;$j<$ntaxa;$j++){
			$parts{"t$i"}{$taxa[$j]} = $bl ? $term_bls[$j] : 0; # stores each terminal branch partitions (id: taxon name)
		}

	}
	my $ntaxa = scalar(keys(%{$parts{t1}}));
	die("Trees '$_[0]' and '$_[1]' are not comparable: different sets of taxa.\n") if $ntaxa != scalar(keys(%{$parts{t2}}));

	my $PD = 0;
	my @t = ("t1","t2");
	for(my $i=0; $i<2; $i++){
		foreach(keys(%{$parts{$t[$i]}})){
			if(exists($parts{$t[$i-1]}{$_})){ # checks if partitions are shared by both trees
				do {$PD += ($parts{$t[$i]}{$_}-$parts{$t[$i-1]}{$_})**2;$parts{$t[$i]}{$_} = 0;$parts{$t[$i-1]}{$_} = 0} if $lensense;
			}
			else{
				die("Trees '$_[0]' and '$_[1]' are not comparable: different sets of taxa.\n") if $_ !~ /,/;
				$PD += $lensense ? $parts{$t[$i]}{$_}**2 : 1;
			}
		}
	}
	$lensense ? return(sqrt($PD)) : $PD;
}
# ============================================================ #


# ====================== AU SUBROUTINES ====================== #
#			(all functions adapted from CONSEL program)
sub getAUpv{
	my $tree = shift; # tree id
	my @rr = @{shift(@_)}; # vector of rescaling rates (== scale factor)
	my @X = @{shift(@_)}; # matrix of predictors for least squares fit
	my @npv = @{shift(@_)}; # vector of naive pvalues (BP)
	my $bb = pop; # n of replicates

	# rcalpcal #	
	my $i;
	my $kk = scalar(@rr); # n of scale factors
	for($i=0;$i<$kk;$i++){
		last if $rr[$i] >= 0;
	}
	my $i1 = $i;
	for(;$i<$kk;$i++){
		last if $rr[$i] >= 1;
	}
	my $i2 = $i;
	
	my $deg; # degeneration switch (j in consel.c)
	my $kappa = 1; # curvature weight
	my ($D,$C,$rss,$se,$pv);
	my ($df,@zval,@wt,@vmat,$pf); # degrees of freedom, vector of response, vector of weights (for weighted least squares) and... pfit(?)

	for($i=$i1;$i<$i2;$i++){
		# wlscalcpval #
		my $j;
		my $m = my $m2 = 0;
		for($j=0;$j<$kk;$j++){
			if($rr[$j]<$rr[$i] || $rr[$j]>100){
				$zval[$j] = $wt[$j] = 0;
			}
			elsif($npv[$j]>=1 || $npv[$j]<=0){
				$zval[$j] = $wt[$j] = 0;
				$m2++;
			}
			else{
				$zval[$j] = -inv_cdf($npv[$j]);
				$wt[$j] = $bb*pdf($zval[$j])**2/((1-$npv[$j])*$npv[$j]);
				$m++;
			}
		}
		$df = $m-2;
		my $x = 0;
		if($m>=2){
			my @fit = lsfit(\@X,\@zval,\@wt,\@vmat); # returns vector containing $D, $C and $rss
			$D = shift(@fit);
			$C = shift(@fit);
			$rss = shift(@fit);
			$pv = cdf(-($D-$kappa*$C));
			$x = pdf($D-$kappa*$C);
			$se = sqrt($x**2*($vmat[0][0]+$kappa**2*$vmat[1][1]-$kappa*$vmat[0][1]-$kappa*$vmat[1][0]));
			$deg = 0;
		}
		else{
			for($j=0;$j<$kk;$j++){
				$x += $zval[$j];
			}
			$pv = ($x>=0) ? 0 : 1;
			$D = $C = $rss = $se = 0;
			$deg = 1;
		}

		##########
		last if $deg;
		$pf = pochisq($rss,$df);
		last if $pf>=0;
	}
	#########

	if($deg){
		my $pf = 0;
		warn("Regression degenerated for $tree: df=$df\n"); # df = negative integer
		$D = $C = 0;
	}
	else{
		if($pf < 0.01){
			warn("Theory does not fit well for $tree: pfit=$pf\n");
		}
	}

	return($pv);
}

sub lsfit{
	my @X = @{shift(@_)}; # matrix of predictors (m x n)
	my @Y = @{shift(@_)}; # vector of response
	my @W = @{shift(@_)}; # vector of weights
	my $vmat = pop;

	my $m = scalar(@X);
	my $n = scalar(@{$X[0]});

	my ($i,$j,$k,$x);

	my (@covmat,@xyvec,@invmat);
	$covmat[0] = [(0,0)]; $covmat[1] = [(0,0)];
	for($i=0;$i<$m;$i++){
		for($x=0,$k=0;$k<$n;$k++){
			$x += $X[$i][$k]*$Y[$k]*$W[$k];
		}
		$xyvec[$i] = $x;
		for($j=0;$j<=$i;$j++){
			for($x=0,$k=0;$k<$n;$k++){
				$x += $X[$i][$k]*$X[$j][$k]*$W[$k];
			}
			$covmat[$j][$i] = $covmat[$i][$j] = $x;
		}
	}

	#printmat(\@covmat);

	@invmat = luinverse(\@covmat);
	#printmat(\@invmat);
	my $sad = sym_mat(\@invmat);
	#printmat(\@invmat);
	warn "lsfit: covariance matrix singularity is $sad" if $sad > 1e-5;
	warn "lsfit: COVARIANCE MATRIX IS SINGULAR" if $sad > 1e-3;
	@{$vmat} = @invmat;

	my @beta = (0,0); # vector for D and C
	for($i=0;$i<$m;$i++){
		for($x=0,$j=0;$j<$m;$j++){
			$x += $invmat[$i][$j]*$xyvec[$j];
		}
		$beta[$i] = $x;
	}

	my $rss = 0;
	my $z;
	for($x=0,$k=0;$k<$n;$k++){
		$z = $X[0][$k]*$beta[0] + $X[1][$k]*$beta[1];
		$rss += $W[$k]*($Y[$k]-$z)**2;
	}

	return (@beta,$rss);
}

sub luinverse{
	my $orimat = shift(@_);
	my $size = scalar(@$orimat);

	my ($i,$j,$k,$maxb);
	my @wk;

	for($i=0; $i<$size; $i++){
		$maxb = 0;
		for($j=0; $j<$size; $j++){
			$maxb = abs($$orimat[$i][$j]) if abs($$orimat[$i][$j])>$maxb;
		}
		die "luinverse: singular matrix" if $maxb == 0;
		$wk[$i] = 1/$maxb;
	}

	my @index;
	my $sum;
	for($j=0; $j<$size; $j++){
		for($i=0; $i<$j; $i++){
			$sum = $$orimat[$i][$j];
			for($k=0; $k<$i; $k++){
				$sum -= $$orimat[$i][$k]*$$orimat[$k][$j];
			}
			$$orimat[$i][$j] = $sum;
		}
		$maxb = 0;
		my $maxi = 0;
		my $tmp;
		for($i=$j; $i<$size; $i++){
			$sum = $$orimat[$i][$j];
			for ($k=0; $k<$j; $k++){
				$sum -= $$orimat[$i][$k]*$$orimat[$k][$j];	
			}
			$$orimat[$i][$j] = $sum;
			$tmp = $wk[$i]*abs($sum);
			if ($tmp >= $maxb){
				$maxb = $tmp;
				$maxi = $i;
			}
		}
		if($j != $maxi){
			for ($k=0; $k<$size; $k++){
				$tmp = $$orimat[$maxi][$k];
				$$orimat[$maxi][$k] = $$orimat[$j][$k];
				$$orimat[$j][$k] = $tmp;
			}
			$wk[$maxi] = $wk[$j];
		}
		$index[$j] = $maxi;
		$$orimat[$j][$j] = 1e-20 if $$orimat[$j][$j] == 0;
		if($j != $size-1){
			$tmp = 1/$$orimat[$j][$j];
			for($i=$j+1; $i<$size; $i++){
				$$orimat[$i][$j] *= $tmp;
			}
		}
	}

	my @invmat = ([()],[()]);
	for(my $jx=0; $jx<$size; $jx++){
		my $ix;
		for($ix=0; $ix<$size; $ix++){
			$wk[$ix] = 0;
		}
		$wk[$jx] = 1.0;
		my $l = -1;
		for($i=0; $i<$size; $i++){
			my $idx = $index[$i];
			$sum = $wk[$idx];
			$wk[$idx] = $wk[$i];
			if($l != -1){
				for($j=$l; $j<$i; $j++){
					$sum -= $$orimat[$i][$j]*$wk[$j];
				}
			} 
			elsif($sum != 0){
				$l = $i;
			}
			$wk[$i] = $sum;
		}
		for($i=$size-1; $i>=0; $i--){
			$sum = $wk[$i];
			for ($j=$i+1; $j<$size; $j++){
				$sum -= $$orimat[$i][$j]*$wk[$j];
			}
			$wk[$i] = $sum/$$orimat[$i][$i];
		}
		for($ix=0; $ix<$size; $ix++){
			$invmat[$ix][$jx] = $wk[$ix];
		}
	}

	return(@invmat);
}

sub sym_mat{
	my $mat = shift(@_);
	my $m = scalar(@{$mat});

	my ($i,$j);
	my $sum = 0;
	for($i=0; $i<$m; $i++){
		for($j=0; $j<$i; $j++){
			$sum += abs($$mat[$i][$j]-$$mat[$j][$i]);
			$$mat[$i][$j] = $$mat[$j][$i] = ($$mat[$i][$j]+$$mat[$j][$i])/2;
		}
	}

	return $sum;
}

sub pochisq{
	my $x = shift; # chi-square value
	my $df = shift; # degrees of freedom

	my $LOG_SQRT_PI = 0.5723649429247000870717135;
	my $I_SQRT_PI = 0.5641895835477562869480795;
	
	return 1 if($x<=0 || $df<1);
	my $a = 0.5*$x;
	my $even = ($df%2 == 0) ? 1 : 0;
	my $y = cexp(-$a) if $df > 1;
	my $s = ($even ? $y : (2*ccdf(-sqrt($x))));

	if($df>2){
		my $x = 0.5*($df-1);
		my $z = ($even ? 1 : 0.5);
		if($a > 20){
			my $e = $even ? 0 : $LOG_SQRT_PI;
			my $c = log($a); 
			for(;$z<=$x;$z++){
				$e = log($z)+$e;
				$s += cexp($c*$z-$a-$e);
			}
			return $s;
		}
		else{
			my $e = $even ? 1 : ($I_SQRT_PI/sqrt($a));
			my $c = 0;
			for(;$z<=$x;$z++){
				$e = $e*($a/$z);
				$c = $c+$e;
			}
			return ($c*$y+$s);
		}
	}
	else{
		return $s;
	}

}

sub cexp{ # conditional exp: returns 0 if x < -20 
	my $x = shift;
	($x < -20) ? return 0 : return exp($x);
}

sub ccdf{ # conditional cdf: returns 0 if abs(x) > 6 
	my $x = shift;
	(abs($x) > 6) ? return 0 : return cdf($x);
}

# ============================================================ #

# =========================== MISC =========================== #

sub removeExt{
	my $in = shift;

	my @temp = split(/\./,$in);
	my $templ = scalar(@temp);

	$templ == 1 ? return $in : return join(".",@temp[0..$templ-2]);
}

sub setSeqAlign{
	my ($hashref,$file_name) = @_;

	my $format;
	my ($nseqs,$nsites) = 0;

	open IN, "<$file_name" or die "Failed to open sequence alignment '$file_name': $!";
	my $header;
	my @seq;

	my $line1 = <IN>;
	until($line1 !~ /^(\h+)?\n/){
		$line1 = <IN>;
	}
	if($line1 =~ /\d+\h+\d+/){ #se for phylip...
		$format = "Phylip";
		while(<IN>){
			if(/(\w+)\h+([\w\?\*-]+)/){
				$header = $1;
				@seq = split("",$2);
				$hashref->{$header} = [@seq];
				$nseqs++;
			}
		}
	}
	elsif($line1 =~ /^>(\w+)/){ #se for fasta...
		$format = "Fasta";
		$header = $1;
		$nseqs++;
		while(<IN>){
			if(/(^[A-Z\?\*-]+)/i){
				push(@seq,split("",$1));
			}
			elsif(/^>(\w+)/){
				$hashref->{$header} = [@seq];
				$header = $1;
				undef @seq;
				$nseqs++;
			}
		}
		$hashref->{$header} = [@seq];
	}
	else{ #se não for nenhum dos dois...
		die "Failed to recognize sequence alignment.\nPlease, make sure it is in supported format (either fasta or sequential phylip).\n";
	}
	close IN;

	$nsites = scalar(@seq);

	return ($nseqs,$nsites,$format);
}

sub setSubstParams{
	my ($hashref,$file_name) = @_;

	#armazena parâmetros de susbtituição estimados
	my @rates;
	my @freqs;
	$hashref->{pinv} = "";
	$hashref->{alpha} = "";
	$hashref->{spr} = "";
	open IN, "<$file_name.iqtree" or die "Failed to open '$file_name.iqtree': $!";
	while(<IN>){
		chomp($_);
		if(/[ACTG]-[ACTG]: (\d\.\d+)/){
			push(@rates,$1);
		}
		elsif(/equal frequencies/){
			@freqs = (0.25) x 4;
		}
		elsif(/pi\([ACTG]\) = (\d\.\d+)/){
			push(@freqs,$1);
		}
		elsif(/Proportion of invariable sites: (\d\.\d+)/){
			$hashref->{pinv} = $1;
		}
		elsif(/Gamma shape alpha: (\d\.\d+)/){
			$hashref->{alpha} = $1;
			last;
		}
		elsif(/Site proportion and rates:\h+(.+)/){
			my $temp = $1;
			$temp =~ s/[\(\)]//g;
			$temp =~ s/,/ /g;
			my @spr = split(" ",$temp);
			$hashref->{spr} = join(",",multiround([@spr],3));
			last;
		}
	}
	close IN;

	$hashref->{freq} = join(",",@freqs);
	#armazena rates no formato utilizado no iqtree
	$hashref->{iqrate} = join(",",uniq(grep($_!=1,@rates)));
	#armazena rates no formato para restrição do modelo GTR no seq-gen
	$hashref->{sgrate} = join(" ",@rates);
	#print defined($hashref->{sgrate});

}

sub setTree{
	my ($hashref,$file_name) = @_;

	open IN, "<$file_name.treefile" or die "Failed to open '$file_name.treefile': $!";
	$hashref->{tree} = <IN>;
	$hashref->{tree} = unroot($hashref->{tree});
	close IN;
}

sub setLikelihoods{
	my ($hashref,$file_name) = @_;

	open IN, "<$file_name.sitelh" or die "Failed to open '$file_name.sitelh': $!";
	seek(IN,length(<IN>),0);
	if(<IN> =~ /^Site_Lh\h+(.+)/){
		$hashref->{slik} = [split(/\h/,$1)];
		$hashref->{lik} = round(sum(@{$hashref->{slik}}),4);
	}
	close IN;
}

sub setLikelihood{
	my ($hashref,$file_name) = @_;

	open IN2, "<$file_name.sitelh" or die "Failed to open '$file_name.sitelh': $!";
	seek(IN2,length(<IN2>),0);
	if(<IN2> =~ /^Site_Lh\h+(.+)/){
		$hashref->{lik} = round(sum(split(/\h/,$1)),4);
	}
	close IN2;
}

sub allEqual{
	my $arref = shift;

	scalar(uniq(@{$arref}))==1 ? return 1 : return 0;
}

sub clean{ # removes all files with pattern given as first argument, unless an exception is specified as second argument
	my $patrn = shift;
	my $excpt = shift if @_;

	if($excpt){
		system('for i in ${patrn}; do if [[ $i != ${excpt} ]]; then rm $i; fi; done');
	}
	else{
		system('rm ${patrn}');
	}
}

sub tnsort{ # uses Schwartzian transformation to emulate a natural sort for tree ids

	return map $_->[0], sort {$a->[2] <=> $b->[2]} map [$_, split /t/], @_;
}

# ============================================================ #

# ================== BEGIN OF MAIN PROGRAM =================== #

#recebe argumentos e identifica, pelos flags, os nomes de arquivos e configurações necessárias aos testes
my ($seq_file,$subst_model,$tree_file,$nreps,@testconf,$seed);
my $verb = 0;
my $ncores = 2;
my $redo = "-redo";
print "\nReading input...\n";
for(my $i=0; $i<scalar(@ARGV); $i+=2){
	if($ARGV[$i] eq "-s"){
		$seq_file = $ARGV[$i+1];
		print "\tSequence alignment file set to '$seq_file'\n";
	}
	elsif($ARGV[$i] eq "-m"){
		$subst_model = $ARGV[$i+1];
		print "\tSubstitution model set to '$subst_model'\n";
	}
	elsif($ARGV[$i] eq "-z"){
		$tree_file = $ARGV[$i+1];
		print "\tTrees file set to '$tree_file'\n";
	}
	elsif($ARGV[$i] eq "-n"){
		$nreps = $ARGV[$i+1];
		print "\tNumber of replicates set to '$nreps'\n";
	}
	elsif($ARGV[$i] eq "-t"){
		@testconf = split("/",$ARGV[$i+1]);
		print "\tTests set to '$ARGV[$i+1]'\n";
	}
	elsif($ARGV[$i] eq "-g"){
		$seed = $ARGV[$i+1];
		print "\tUser provided seed: $ARGV[$i+1]\n";
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
			warn "\tWarning: Reusing any previously computed likelihood!!!\n";
		}
		elsif($ARGV[$i+1] !~ /^T(RUE)?/){
			die "Invalid entry '$ARGV[$i+1]' for redoing mode (-redo). It must be a boolean (TRUE/FALSE, or shortly T/F). Please, reset and try again.\n";
		}
	}
	else{
		die "Unrecognized flag '$ARGV[$i]': Check for possible typos or missing values.\n\nAll allowed flags and corresponding values:\n\t-s <sequence_alignment_file_name> (fasta or phylip)\n\t-m <substitution_model> (iqtree format)\n\t-t <trees_file> (newick format - one per line)\n\t-n <number_of_replicates> (default = 1000)\n\t-t <tests_demanded> (default = BP/KH/SH)\n\nFor more info, please check the manual.\n";
	}
}
if(!defined($seq_file)||!defined($subst_model)||!defined($tree_file)){
	print "Mandatory argument(s) missing:\n";
	if(!defined($seq_file)){
		print "\t(-s) sequence alignment file name\n";
	}
	if(!defined($subst_model)){
		print "\t(-m) substitution model\n";
	}
	if(!defined($tree_file)){
		print "\t(-z) trees file name\n";
	}
	die;
}
if(!defined($nreps)){
	print "\tNumber of replicates (-n) not specified. Setting to default '1000'\n";
	$nreps = 1000;
}
if(!@testconf){
	print "\tTopology tests (-t) not specified. Setting to default 'BP/KH/SH'\n";
	@testconf = ("BP","KH","SH");
}
if(!$seed){
	$seed = time();
	print "\tTime based seed: $seed\n";
}



#armazena os testes requisitados e os procedimentos pelos quais devem ser realizados (aproximação normal, bootstrap paramétrico ou não-paramétrico e níveis de otimização dos parâmetros)
my ($BP,$KH,$SOWH,$SH,$ELW,$AU) = 0; # test switches
my %pt; # hash of 'Procedure vs Test' switches
foreach(@testconf){
	my @temp = split(/:/);
	if($temp[0] eq "BP"){
		$BP = 1;
	}
	elsif($temp[0] eq "KH"){
		$KH = 1;
	}
	elsif($temp[0] eq "SOWH"){
		$SOWH = 1;
	}
	elsif($temp[0] eq "SH"){
		$SH = 1;
	}
	elsif($temp[0] eq "ELW"){
		$ELW = 1;
	}
	elsif($temp[0] eq "AU"){
		$AU = 1;
	}
	else{
		die "'$temp[0]' is an invalid test option.\nPlease, try one of the tests available in phyTest (BP, AU, KH, SH, ELW or SOWH).\n";
	}
	if(defined($temp[1])){
		my @p = sort(split(/,/,$temp[1]));
		foreach(@p){
			if($_<-3 || $_>4){
				die "'$_' is an invalid procedure option.\nPlease, try one of the procedures available in phyTest (-3,-2,-1,0,1,2,3,4). Consult manual for further details.\n";
			}
			if(($_ == -3 && $temp[0] ne "KH") || ($_>-3 && $_<0 && $temp[0] ne "SOWH") || ($_>=0 && $_ <=4 && $temp[0] eq "SOWH") || ($_>2 && ($temp[0] ne "AU" && $temp[0] ne "BP"))){
				die "Procedure parameter '$_' is not compatible with $temp[0] test.\nPlease, consult manual for a valid procedure, reset required tests argument (-t) and try again.\n";
			}
			$pt{$_}{$temp[0]} = 1;
		}
	}
	elsif($temp[0] eq "SOWH"){
		$pt{-1}{SOWH} = 1;
	}
	else{
		$pt{0}{$temp[0]} = 1;
	}
}
my $pnum = scalar(keys(%pt)); # number of procedures

my @p; # array of procedure switches
$p[7] = "";
foreach(keys(%pt)){
	$p[$_] = 1;
}



#armazena o modelo de substituição a ser utilizado para calcular todas as verossimilhanças
if($subst_model=~/\{.+\}/){
	die "This version of phyTest cannot handle pre-specified substitution parameter values (any between '{}') in '$subst_model'.\nPlease, reset the model argument (-m) and try again.\n";
}
my @params = split(/\+/,$subst_model);
my $matrix = shift @params;
my $inv = "";
my $gama = "";
my $free = "";
foreach(@params){
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
	die "Cannot use '+G' and '+R' simultaniously as passed in '$subst_model'.\nPlease, reset the model argument (-m) and try again.\n";
}
elsif($free && $SOWH){
	die "Due to limitations in Seq-Gen program, freely distributed rate categories (+R) - as passed in '$subst_model' - can not be used with SOWH test option.\nPlease, reset the model (-m) or required tests (-t) and try again.\n";
}
undef @params;



my $data_name = removeExt($seq_file);
my $nsites;
my %data;
print "\nReading sequence alignment from '$seq_file'...\n";
#se for necessário reamostrar sitios do alinhamento (dado em $seq_file) para bootstrap não-paramétrico, identifica seu formato (phylip ou fasta) e o armazena em %data, além do seu número de sítios em $nsites
if($p[1]||$p[2]||$p[3]||$p[4]){

	my @info = setSeqAlign(\%data,$seq_file);

	print "\t$info[2] format detected\n";
	print "\t$info[0] sequences\n";
	$nsites = $info[1];
	print "\t$nsites sites long\n";

}
# caso contrário, se for realizar bootstrap paramétrico ou RELL, obtém apenas o número de sítios pelo tamanho da primeira sequencia
elsif($p[0]||$p[-1]||$p[-2]){

	open IN, "<$seq_file";
	while(<IN>){
		if(/^(\w+\h+)?([\w\?\*-]+)/){
			if(defined($1)){
				print "\tPhylip format detected\n";
			}
			else{
				print "\tFasta format detected\n";
			}
			$nsites = length($2);
			print "\t$nsites sites long\n";
			last;
		}
	}
	close IN;
}



#recebe as árvores a serem testadas (em formato newick e listadas em um mesmo doc de texto), associa cada uma a uma ID e as desenraiza
my @trees;
open IN, "<$tree_file" or die "Failed to open '$tree_file': $!";
print "\nReading phylogenetic trees from '$tree_file'...\n";
my $id;
while(<IN>){
	if(/[\(\)\w\.,:]+;/){
		$id++;
		push(@trees,unroot($_));
	}
}
close IN;
# garante que as árvores recebidas sejam comparáveis, topologicamente únicas e que o número de topologias únicas seja suficiente para os testes requisitados; armazena as árvores em $opt{d}, gerando um arquivo separado para cada uma
my %opt;
my $ntrees = scalar(@trees); #total number of trees
my $nutrees; # number of topologically unique trees
if($ntrees==0){
	die "\nNot enough trees to proceed: no phylogenetic tree detected.\n";
}
else{
	print "\t$ntrees tree(s) detected and successfully unrooted. ";
	if($ntrees==1){
		open OUT, ">t1.tre";
		print OUT $trees[0];
		close OUT;
		$opt{d}{t1}{tree} = $trees[0];
	}
	elsif($ntrees>1){
		my ($i,$j);
		for($i=0; $i<$ntrees; $i++){
			if(defined($trees[$i])){
				$id = "t".($i+1);
				open OUT, ">${id}.tre";
				print OUT $trees[$i];
				close OUT;
				$opt{d}{$id}{tree} = $trees[$i];
				for($j=$i+1; $j<$ntrees; $j++){
					if(defined($trees[$j]) && phyDist($trees[$i],$trees[$j])==0){ # ao buscar diferenças topológicas entre as árvores, a função phyDist também garante que sejam comparáveis (mesmo set de táxons)
						$trees[$j] = undef;
					}
				}
			}
		}
		undef @trees;
		@trees = tnsort(keys(%{$opt{d}})); # array of topologically unique trees
		$nutrees = scalar(@trees);

		if($ntrees==$nutrees){
			print "All topologically unique.\n"
		}
		else{
			print "Some have the same topology.\n\tOnly the following $nutrees will be considered: ".join(", ",@trees)."\n";
		}
	}
	if($ntrees==1 || $nutrees==1){
	 	if($pnum>2 || ($pnum==2 && !$p[3] && !$p[4]) || ($pnum==1 && !$p[3] || !$p[4])){
			die "\nNot enough trees to proceed: using the tests and/or procedures passed via (-t), at least 2 reasonable and topologically unique trees must be provided to achieve minimal precision in confidence calculation.\n";
		}
	}	
}



#dado o alinhamento original, utiliza o iqtree para otimizar os comprimentos de ramo, além dos parâmetros de substituição para cada árvore e os armazena junto a suas verossimilhanças em %opt
my $dmlt;
foreach my $t(@trees){ # para cada árvore $t...
	print "\n\t$t ";
	#via iqtree, otimiza seus comprimentos de ramo e parametros de substituição para o alinhamento dado
	system "iqtree-omp -s $seq_file -te ${t}.tre -m $subst_model -pre ${t}_${data_name}_${matrix}${gama}${free}${inv} -nt $ncores -wsl -quiet -redo" if (!(-f "${t}_${data_name}_${matrix}${gama}${free}${inv}.sitelh") || !(-f "${t}_${data_name}_${matrix}${gama}${free}${inv}.treefile") || !(-f "${t}_${data_name}_${matrix}${gama}${free}${inv}.iqtree") || $redo);

	#rearmazena árvore com comprimentos de ramo otimizados
	setTree(\%{$opt{d}{$t}},"${t}_${data_name}_${matrix}${gama}${free}${inv}");

	#armazena valor de log-verossimilhança total da árvore e suas log-verossimilhanças por sítio
	setLikelihoods(\%{$opt{d}{$t}},"${t}_${data_name}_${matrix}${gama}${free}${inv}");

	print "lnL = $opt{d}{$t}{lik}";
	# identifica árvore de maior log-verossimilhança para o alinhamento ($dmlt), entre as árvores dadas
	if(!defined($dmlt) || $opt{d}{$dmlt}{lik} < $opt{d}{$t}{lik}){
		$dmlt = $t;
	}

}
print "\n\n\tHighest-scoring tree: $dmlt (when performing KH, SH or SOWH test, the remaining trees are compared against this)\n";
if($KH||$SH||$SOWH){ # se for realizado qualquer teste baseado em delta de verossimilhança...
	# armazena diferenças de log-verossimilhaca entre a $dmlt e demais árvores
	foreach my $t(@trees){
		$opt{d}{$t}{delta} = round($opt{d}{$dmlt}{lik}-$opt{d}{$t}{lik},4);
	}
}

if($p[-1]||$p[-2]||$p[1]||$p[2]){ # se for realizar qualquer bootstrap sem busca por MLT...

	foreach my $t(@trees){
		setSubstParams(\%{$opt{d}{$t}},"${t}_${data_name}_${matrix}${gama}${free}${inv}");
	}

}
my %max;
if($p[3] || $p[4]){ # se for necessária busca pelos parâmetros de subs. de máxima verossimilhança...

	print "\nSearching for maximum likelihood tree for ${data_name}...";
	system "iqtree-omp -s $seq_file -m $subst_model -pre MLt_${data_name}_${matrix}${gama}${free}${inv} -nt $ncores -wsl -quiet -redo" if (!(-f "MLt_${data_name}_${matrix}${gama}${free}${inv}.sitelh") || $redo);

	setSubstParams(\%{$max{d}},"MLt_${data_name}_${matrix}${gama}${free}${inv}");

}

my $Coptfreq;
if($p[-2]||$p[2]||$p[4]){ # se for realizar qualquer tipo de bootstrap com otimização completa...
	my $freqs = $p[4] ? $max{d}{freq} : $opt{d}{$dmlt}{freq};

	if(allEqual([split(",",$freqs)])){
		$Coptfreq = "+FQ";
	}
	else{
		$Coptfreq = "+FO";
	}
}



my %pvs;
if($p[0]||$p[1]||$p[2]||$p[3]||$p[4]){ # se o(s) teste(s) envolver(em) geração de réplicas não-paramétricas...

	srand($seed);

	my @procs = grep($_>=0,keys(%pt)); # lista de procedimentos
	# define uma lista de fatores de escala, para caso venha a ser realizado o teste AU, ou apenas um fator redundante (1), caso não
	my @sf = $AU ? (0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4) : (1);


	for(my $k=0; $k<scalar(@sf); $k++){ # para cada fator de escala...

		my $snsites = round($nsites*$sf[$k],0); # scaled number of sites
		print "\n\nGenerating $nreps non-parametric replicates with $snsites sites each and executing any required optimization...\n\n";
		print "  0%................................................100%\n";
		print "   ";
		for(my $r=1; $r<=$nreps; $r++){ # dada cada réplica $r...

			#sorteia, com reposição, as posições a serem amostradas para compor a réplica $r, totalizando N' posições = N posições do alinhamento original ($nsites) X o fator de escala ($sf[$k])
			my @posits;
			for(my $i=0; $i<$snsites; $i++){
				$posits[$i] = int(rand($nsites));
			}

			#gera alinhamento-réplica, se necessário para algum dos procedimentos
			if($p[1]||$p[2]||$p[3]||$p[4]){

				my $rep;
				foreach my $header(keys(%data)){
					$rep .= ">".$header."\n".join("",@{$data{$header}}[@posits])."\n";
				}
				open(OUT,">","${data_name}_NPboot${r}.fas");
				print OUT $rep;
				close OUT;

			}

			foreach(@procs){ # para cada procedimento $_...

				if($sf[$k] == 1 || $pt{$_}{AU}){ # se o fator de escala atual for cabível,
					my $rmlt;

					if($_ < 3){ # e se o procedimento não envolver busca pela árvore ML,					
						foreach my $t(@trees){ # para cada árvore $t...

							#RELL
							if($_ == 0){
								print "\tResampling estimated log-likelihoods (RELL) for each tree: replicate $r (".$nsites*$sf[$k]." sites long)\n" if $verb && $t eq "t1";
								#reamostra os valores de verossimilhanças por sítio para cada árvore e armazena os totais resultantes
								$opt{r}{$t}{lik} = round(sum(@{$opt{d}{$t}{slik}}[@posits]),4);

							}
							#Bootstrap não-paramétrico
							else{

								# utiliza o alinhamento-réplica para otimizar os parâmetros de cada árvore
								if($_ == 1){ # com otimização parcial (apenas comprimentos de ramos)...
									print "\tOptimizing trees using replicate $r (".$nsites*$sf[$k]." sites long)\n" if $verb && $t eq "t1";
									system "iqtree-omp -s ${data_name}_NPboot${r}.fas -te ${t}.tre -m $matrix'{'$opt{d}{$t}{iqrate}'}'$inv'{'$opt{d}{$t}{pinv}'}'$gama'{'$opt{d}{$t}{alpha}'}'$free'{'$opt{d}{$t}{spr}'}'+F'{'$opt{d}{$t}{freq}'}' -pre ${t}_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_Popt -nt $ncores -wsl -quiet -redo" if (!(-f "${t}_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_Popt.sitelh") || $redo);
									
									setLikelihood(\%{$opt{r}{$t}},"${t}_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_Popt");
								}
								else{ # ou completa (ramos + parametros de susbtituição)
									print "\tOptimizing trees and substitution parameters using replicate $r (".$nsites*$sf[$k]." sites long)\n" if $verb && $t eq "t1";
									system "iqtree-omp -s ${data_name}_NPboot${r}.fas -te ${t}.tre -m ${matrix}${inv}${gama}${free}${Coptfreq} -pre ${t}_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_Copt -nt $ncores -wsl -quiet -redo" if (!(-f "${t}_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_Copt.sitelh") || $redo);

									setLikelihood(\%{$opt{r}{$t}},"${t}_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_Copt");
								}
								#clean("${t}_${data_name}_NPboot${r}*","*.sitelh");

							}
							# identifica a árvore de maior verossimilhança para a réplica
							if(!defined($rmlt) || $opt{r}{$rmlt}{lik} < $opt{r}{$t}{lik}){
								$rmlt = $t;
							}

						}
					}
					else{ # se o procedimento envolver busca pela árvore ML,
						# utiliza o alinhamento-réplica para inferí-la
						if($_ == 3){ # com otimização parcial (apenas topologia + comprimentos de ramos)
							print "\tSearching ML tree for replicate $r (".$nsites*$sf[$k]." sites long)\n" if $verb;
							system "iqtree-omp -s ${data_name}_NPboot${r}.fas -m $matrix'{'$max{d}{iqrate}'}'$inv'{'$max{d}{pinv}'}'$gama'{'$max{d}{alpha}'}'$free'{'$max{d}{spr}'}'+F'{'$max{d}{freq}'}' -pre MLt_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_TrSch_Popt -nt $ncores -wsl -quiet -redo" if (!(-f "MLt_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_TrSch_Popt.treefile") || $redo);
							
							open IN, "<MLt_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_TrSch_Popt.treefile" or die "Failed to open 'MLt_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_TrSch_Popt.treefile': $!";
							$rmlt = <IN>;
							close IN;
						}
						else{ # ou completa (topologia + ramos + parametros de substituição)
							print "\tSearching ML tree and substitution parameters for replicate $r (".$nsites*$sf[$k]." sites long)\n" if $verb;
							system "iqtree-omp -s ${data_name}_NPboot${r}.fas -m ${matrix}${inv}${gama}${free}${Coptfreq} -pre MLt_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_TrSch_Copt -nt $ncores -wsl -quiet -redo" if (!(-f "MLt_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_TrSch_Copt.treefile") || $redo);
							
							open IN, "<MLt_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_TrSch_Copt.treefile" or die "Failed to open 'MLt_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_TrSch_Copt.treefile': $!";
							$rmlt = <IN>;
							close IN;
						}
						$rmlt = unroot($rmlt);
						#clean("MLt_${data_name}_${matrix}${gama}${free}${inv}_NPboot${r}_TrSch*","*.treefile");
					}

					my $elw_ratiosum = 0 if $sf[$k] == 1 && $pt{$_}{ELW};
					foreach my $t(@trees){

						if($sf[$k] == 1){
							if($pt{$_}{KH}){ # se foi requisitado KH sob o procedimento
								push(@{$pvs{$t}{KH}{$_}{rdeltas}},$opt{r}{$dmlt}{lik}-$opt{r}{$t}{lik}); # armazena delta entre a árvore que foi a de ML para o alinhamento original e $t (para posterior centralização)
							}
							if($pt{$_}{SH}){ # se foi requisitado SH
								push(@{$pvs{$t}{SH}{$_}{rliks}},$opt{r}{$t}{lik}); # armazena apenas lnL de $t (para posterior centralização)
							}
							if($pt{$_}{ELW}){ # se foi requisitado ELW
								$pvs{$t}{ELW}{$_}{rratio} = exp($opt{r}{$t}{lik}-$opt{r}{$rmlt}{lik}); # armazena razão entre a lnL de $t e da árvore de ML para $r (para computação mais precisa da razão, é obtida de e^delta, sendo delta = diferença entre lnLs)
								$elw_ratiosum += $pvs{$t}{ELW}{$_}{rratio}; # e a acumula à soma das razões para as diferentes árvores
							}
						}
						if($pt{$_}{BP}||$pt{$_}{AU}){ # se foi requisitado BP ou AU sob o procedimento,
							$pvs{$t}{BP}{$_}{$k} = 0 if !exists($pvs{$t}{BP}{$_}{$k});
							if($_ < 3 && $t eq $rmlt){ # sendo ele RELL, otimização parcial ou total e se $t tiver obtido maior lnL para a réplica atual,
								$pvs{$t}{BP}{$_}{$k}++; # adiciona à contagem de vezes em que $t foi a árvore de ML
							}
							elsif($_ >= 3){ # ou sendo ele busca topológica + otimização parcial ou total e se $t tiver a topologia de ML para a réplica atual,
								if(phyDist($opt{d}{$t}{tree},$rmlt) == 0){
									$pvs{$t}{BP}{$_}{$k}++; # adiciona à contagem de vezes em que $t foi a árvore de ML
								}
							}
						}

					}
					$rmlt = undef;
					# terminadas todas as árvores
					if($sf[$k] == 1 && $pt{$_}{ELW}){ # se foi requisitado ELW sob o procedimento
						foreach my $t(@trees){ # novamente, dada cada árvore $t...
							$pvs{$t}{ELW}{$_}{weightsum} = 0 if !exists($pvs{$t}{ELW}{$_}{weightsum});
							$pvs{$t}{ELW}{$_}{weightsum} += round(round($pvs{$t}{ELW}{$_}{rratio},4)/round($elw_ratiosum,4),4); # armazena peso da razão de sua verossimilhança (ponderado pela soma das razões, calculada acima)
							delete $pvs{$t}{ELW}{$_}{rratio};
						}
						$elw_ratiosum = undef;
					}

				}
			}

			system "rm ${data_name}_NPboot${r}.fas" if $p[1]||$p[2]||$p[3]||$p[4];
			print "|" if $r%round($nreps/50) == 0;

		}
		# terminadas todas as réplicas do fator $k, calcula pvalues (de BP, KH, SH e ELW apenas)
		if($sf[$k] == 1){
			
			foreach(@procs){ # novamente, para cada procedimento $_...

				foreach my $t(@trees){ # e para cada árvore $t...

					if($pt{$_}{BP}){ # se foi requisitado BP

						print "\tComputing confidence on each tree via BP:$_\n" if $verb && $t eq "t1";
						$pvs{$t}{BP}{$_}{pval} = $pvs{$t}{BP}{$_}{$k}/$nreps;

					}
					if($pt{$_}{KH}){ # se foi requisitado KH

						print "\tComputing confidence on each tree via KH:$_\n" if $verb && $t eq "t1";
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
					if($pt{$_}{SH}){ # se for requisitado SH

						$pvs{$t}{SH}{$_}{pval} = 0 if !exists($pvs{$t}{SH}{$_}{pval});
						$pvs{$t}{SH}{$_}{likmean} = round((sum(@{$pvs{$t}{SH}{$_}{rliks}})/$nreps),4); # armazena a média das lnLs entre réplicas temporariamente

					}
					if($pt{$_}{ELW}){ # se foi requisitado ELW

						print "\tComputing confidence on each tree via ELW:$_\n" if $verb && $t eq "t1";
						$pvs{$t}{ELW}{$_}{pval} = round($pvs{$t}{ELW}{$_}{weightsum}/$nreps,4);
						delete $pvs{$t}{ELW}{$_}{weightsum};

					}
				}
				# terminadas todas as árvores
				if($pt{$_}{SH}){ # se foi requisitado SH
					
					for(my $r=1; $r<=$nreps; $r++){ # para cada réplica $r...
						my $rmlt;
						foreach my $t(@trees){ # e para cada árvore $t...
							$pvs{$t}{SH}{$_}{rliks}[$r] -= $pvs{$t}{SH}{$_}{likmean}; # executa o procedimento de centralização de lnLs
							# seleciona a árvore que obtém maior lnL centralizado para a réplica
							if(!defined($rmlt) || $pvs{$rmlt}{SH}{$_}{rliks}[$r] < $pvs{$t}{SH}{$_}{rliks}[$r]){
								$rmlt = $t;
							}
						}
						foreach my $t(@trees){ # para cada árvore, novamente...
							if($opt{d}{$t}{delta}<=($pvs{$rmlt}{SH}{$_}{rliks}[$r]-$pvs{$t}{SH}{$_}{rliks}[$r])){ # se o delta entre e a árvore de ML para $r e $t superar o delta original...
								$pvs{$t}{SH}{$_}{pval}++; # adiciona à contagem de vezes em que o delta original foi superado
							}
						}
					}
					print "\tComputing confidence on each tree via SH:$_\n" if $verb;
					foreach my $t(@trees){
						$pvs{$t}{SH}{$_}{pval} = round($pvs{$t}{SH}{$_}{pval}/$nreps,4);
						delete $pvs{$t}{SH}{$_}{likmean};
						delete $pvs{$t}{SH}{$_}{rliks};
					}

				}

			}

		}

	}
	#terminados todos os fatores de escala, calcula pvalue do AU
	if($AU){

		my $i=0;
		my @X = ([()],[()]);
		foreach(@sf){
			$X[0][$i] = sqrt($sf[$i]); # vector of rooted rescaling rates
			$X[1][$i] = 1/$X[0][$i]; # vector of inverse of rooted rescaling rates
			$i++;
		}

		foreach(@procs){
			print "\tComputing confidence on each tree via AU:$_\n" if $verb;

			foreach my $t(@trees){

				my @npv;
				my $kk = scalar(@sf);
				for(my $k=0;$k<$kk;$k++){
					my $naive_pv = $pvs{$t}{BP}{$_}{$k}/$nreps;
					$npv[$k] = $naive_pv;
					delete($pvs{$t}{BP}{$_}{$k});
				}

				$pvs{$t}{AU}{$_}{pval} = getAUpv($t,\@sf,\@X,\@npv,$nreps);

				delete $pvs{$t}{BP}{$_} if !$pt{$_}{BP};
				delete $pvs{$t}{BP} if !$BP;
			}
		}

	}

}
if($p[-1]||$p[-2]){ # se envolver(em) geração de réplicas paramétricas...

	my @procs = grep(/-(1|2)/,keys(%pt));
	# converte número de categorias gama para o formato do seq-gen
	my $gcat = "";
	if($gama){
		if($gama =~ /\+G(\d+)/){
			$gcat = "-g $1";
		}
		else{
			$gcat = "-g 4";
		}
	}

	foreach my $truet (@trees){ # para cada árvore $truet

		# converte alpha e proporção de sítios invariáveis
		my $alph = "";
		if($gama){
			$alph = "-a $opt{d}{$truet}{alpha}";
		}
		my $pinv = "";
		if($inv){
			$pinv = "-i $opt{d}{$truet}{pinv}";
		}
		# e simula N alinhamentos utilizando a árvore $truet, sendo N = número de réplicas dado em $nreps
		print "\nUsing $truet to simulate $nreps parametric replicates with ".$nsites." sites each\n";
		print "  0%................................................100%\n";
		print "   ";
		system "seq-gen -m GTR -l $nsites -n $nreps $gcat $alph $pinv -f $opt{d}{$truet}{freq} -r $opt{d}{$truet}{sgrate} -or -q -z $seed < ${truet}.tre > ${data_name}-${truet}_Pboots.data";	
		open IN, "<${data_name}-${truet}_Pboots.data" or die "Failed to open '${data_name}-${truet}_Pboots.data': $!";
		local $/ = "\n ";
		my $r;	
		while(<IN> =~ /(\d+\h\d+\n(\w+\h[ACGT]+\n)+)/){ # para cada alinhamento gerado...
			$r++;

			# cria arquivo separado o contendo
			open OUT, ">${data_name}-${truet}_Pboot${r}.phy";
			print OUT $1;
			close OUT;
			foreach(@procs){ # para cada procedimento $_...

				my $rmlt;
				foreach my $t(@trees){ # e, para cada árvore, utiliza o alinhamento simulado para otimização

					if($_ == -1){	# parcial (apenas dos comprimentos de ramos, fixando os mesmos parâmetros de substituição utilizados para a simulação dos dados)
						print "\tOptimizing trees using replicate $r (true tree: $truet)\n" if $verb && $t eq "t1";
						system "iqtree-omp -s ${data_name}-${truet}_Pboot${r}.phy -te ${t}.tre -m $matrix'{'$opt{d}{$truet}{iqrate}'}'$inv'{'$opt{d}{$truet}{pinv}'}'$gama'{'$opt{d}{$truet}{alpha}'}'+F'{'$opt{d}{$truet}{freq}'}' -pre ${t}_${data_name}_${matrix}${gama}${free}${inv}-${truet}_Pboot${r}_Popt -nt $ncores -wsl -quiet -redo" if (!(-f "${t}_${data_name}_${matrix}${gama}${free}${inv}-${truet}_Pboot${r}_Popt.sitelh") || $redo);
						
						setLikelihood(\%{$opt{r}{$t}},"${t}_${data_name}_${matrix}${gama}${free}${inv}-${truet}_Pboot${r}_Popt");
						system "rm ${t}_${data_name}-${truet}_Pboot${r}_Popt*";
					}
					else{	# ou completa (ramos + parametros de susbtituição)
						print "\tOptimizing trees and substitution parameters using replicate $r (true tree: $truet)\n" if $verb && $t eq "t1";
						system "iqtree-omp -s ${data_name}-${truet}_Pboot${r}.phy -te ${t}.tre -m ${matrix}${inv}${gama}${Coptfreq} -pre ${t}_${data_name}_${matrix}${gama}${free}${inv}-${truet}_Pboot${r}_Copt -nt $ncores -wsl -quiet -redo" if (!(-f "${t}_${data_name}_${matrix}${gama}${free}${inv}-${truet}_Pboot${r}_Copt.sitelh") || $redo);
						
						setLikelihood(\%{$opt{r}{$t}},"${t}_${data_name}_${matrix}${gama}${free}${inv}-${truet}_Pboot${r}_Copt");
						system "rm ${t}_${data_name}-${truet}_Pboot${r}_Popt*";
					}
					# e define a árvore de maior verossimilhança
					if(!defined($rmlt) || $opt{r}{$rmlt}{lik} < $opt{r}{$t}{lik}){
						$rmlt = $t;
					}

				}
				$pvs{$truet}{SOWH}{$_}{pval} = 0 if !exists($pvs{$truet}{SOWH}{$_}{pval});
				if($opt{d}{$truet}{delta}<=($opt{r}{$rmlt}{lik}-$opt{r}{$truet}{lik})){ # se o delta entre e a árvore de ML para $r e $truet superar o delta original...
					$pvs{$truet}{SOWH}{$_}{pval}++; # adiciona à contagem de vezes em que o delta original foi superado
				}

			}

			system "rm ${data_name}-${truet}_Pboot${r}.phy";
			print "|" if $r%round($nreps/50) == 0;
		}
		close IN;

		system "rm ${data_name}-${truet}_Pboots.data";

		#terminadas todas as réplicas paramétricas, calcula o pvalue de $truet
		foreach(@procs){ # novamente, para cada procedimento $_...
			print "\tComputing confidence on $truet via SOWH:$_\n" if $verb;
			$pvs{$truet}{SOWH}{$_}{pval} = $pvs{$truet}{SOWH}{$_}{pval}/$nreps;
		}

	}

}
delete $opt{r};
if($p[-3]){ # se envolver(em) presunção de normalidade de delta...

	print "\tComputing confidence on each tree via KH:-3\n";
	foreach my $t(@trees){ # para cada árvore $t...

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
print "\n\n";
printf("%7s%13s","tree_id","lnL");
printf("%10s","delta") if $KH||$SH||$SOWH;
foreach my $test(sort(keys(%{$pvs{t1}}))){
	foreach(sort(keys(%{$pvs{t1}{$test}}))){
		printf("%7s","$test:$_");
	}
}
print "\n";
foreach my $tree(@trees){
	printf("%7s%13.3lf",$tree,$opt{d}{$tree}{lik});
	printf("%10.3lf",$opt{d}{$tree}{delta}) if $KH||$SH||$SOWH;
	foreach my $test(sort(keys(%{$pvs{$tree}}))){
		foreach(sort(keys(%{$pvs{$tree}{$test}}))){
			printf("%7.3f",$pvs{$tree}{$test}{$_}{pval});
		}
	}
	print "\n";
}
print "\n";

# ================= END OF MAIN PROGRAM ================== #

exit;