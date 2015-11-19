drop view if exists tmp_gbsgenotype;
create view tmp_gbsgenotype as(
	select rs, chrom, seqkey.sample, seqkey.plate_well, genotype,
		case when seqkey.plate_well = '3B12' then genotype end as 3B12,
		case when seqkey.plate_well = '3A09' then genotype end as 3A09
	from genotype join seqkey
		on genotype.dna_sample = seqkey.dna_id
		where call_id = 'GBS'
		and seqkey.sample like 'StNY001'
		order by rs
);
drop view if exists tmp_gbsgenotype_pivot;
create view tmp_gbsgenotype_pivot as(
	select 
		rs, chrom,
		group_concat(3B12) as 3B12,
		group_concat(3A09) as 3A09
	from tmp_gbsgenotype
	group by chrom, rs
);
select 'rs', 'chrom', '3B12', '3A09'
UNION ALL
select * from tmp_gbsgenotype_pivot
into outfile '~/workdir/stny001.txt'
;
