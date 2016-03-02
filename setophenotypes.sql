select * from trait;
select trait_id, trait, description from trait;
select trait_id, unit_category, factors, location from trait;

select isolate, subset, block, plot, phenotype
from phenotype
where trait = 'abun3'
order by plot
into outfile '~/workdir/abun3.txt'
;
select isolate, subset, block, plot, phenotype
from phenotype
where trait = 'abun5'
order by plot
into outfile '~/workdir/abun5.txt'
;
select isolate, subset, block, plot, phenotype
from phenotype
where trait = 'ci'
order by plot
into outfile '~/workdir/ci.txt'
;
select isolate, subset, block, dpi, plate, plot, phenotype
from phenotype
where trait = 'diam'
order by plot
into outfile '~/workdir/diam.txt'
;
select isolate, subset, block, genotype, plot, phenotype
from phenotype
where trait = 'dla'
order by plot
into outfile '~/workdir/dla.txt'
;
select isolate, subset, block, genotype, plot, phenotype
from phenotype
where trait = 'dla814'
order by plot
into outfile '~/workdir/dla814.txt'
;
select isolate, subset, block, plot, plant, genotype, phenotype
from phenotype
where trait = 'ip'
order by plot
into outfile '~/workdir/ip.txt'
;
select isolate, subset, block, plot, phenotype
from phenotype
where trait = 'mel'
order by plot
into outfile '~/workdir/mel.txt'
;
select isolate, subset, block, plot, rep, counts, phenotype
from phenotype
where trait = 'sp'
order by plot
into outfile '~/workdir/sp.txt'
;
