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
drop view if exists tmp_gbsgenotype;
create view tmp_gbsgenotype as(
	select rs, chrom, seqkey.sample, seqkey.plate_well, genotype,
		case when seqkey.plate_well = '3B09' then genotype end as 3B09,
		case when seqkey.plate_well = '3C12' then genotype end as 3C12
	from genotype join seqkey
		on genotype.dna_sample = seqkey.dna_id
		where call_id = 'GBS'
		and seqkey.sample like 'St52B'
		order by rs
);
drop view if exists tmp_gbsgenotype_pivot;
create view tmp_gbsgenotype_pivot as(
	select 
		rs, chrom,
		group_concat(3B09) as 3B09,
		group_concat(3C12) as 3C12
	from tmp_gbsgenotype
	group by chrom, rs
);
select 'rs', 'chrom', '3B09', '3C12'
UNION ALL
select * from tmp_gbsgenotype_pivot
into outfile '~/workdir/st52B.txt'
;
drop view if exists tmp_gbsgenotype;
create view tmp_gbsgenotype as(
	select rs, chrom, seqkey.sample, seqkey.plate_well, genotype,
		case when seqkey.plate_well = '3C11' then genotype end as 3C11,
		case when seqkey.plate_well = '3D12' then genotype end as 3D12
	from genotype join seqkey
		on genotype.dna_sample = seqkey.dna_id
		where call_id = 'GBS'
		and seqkey.sample like 'St54A'
		order by rs
);
drop view if exists tmp_gbsgenotype_pivot;
create view tmp_gbsgenotype_pivot as(
	select 
		rs, chrom,
		group_concat(3C11) as 3C11,
		group_concat(3D12) as 3D12
	from tmp_gbsgenotype
	group by chrom, rs
);
select 'rs', 'chrom', '3C11', '3D12'
UNION ALL
select * from tmp_gbsgenotype_pivot
into outfile '~/workdir/st54A.txt'
;
drop view if exists tmp_gbsgenotype;
create view tmp_gbsgenotype as(
	select rs, chrom, seqkey.sample, seqkey.plate_well, genotype,
		case when seqkey.plate_well = '3B10' then genotype end as 3B10,
		case when seqkey.plate_well = '3D06' then genotype end as 3D06
	from genotype join seqkey
		on genotype.dna_sample = seqkey.dna_id
		where call_id = 'GBS'
		and seqkey.sample like 'St34A'
		order by rs
);
drop view if exists tmp_gbsgenotype_pivot;
create view tmp_gbsgenotype_pivot as(
	select 
		rs, chrom,
		group_concat(3B10) as 3B10,
		group_concat(3D06) as 3D06
	from tmp_gbsgenotype
	group by chrom, rs
);
select 'rs', 'chrom', '3B10', '3D06'
UNION ALL
select * from tmp_gbsgenotype_pivot
into outfile '~/workdir/st34A.txt'
;
drop view if exists tmp_gbsgenotype;
create view tmp_gbsgenotype as(
	select rs, chrom, seqkey.sample, seqkey.plate_well, genotype,
		case when seqkey.plate_well = '3C10' then genotype end as 3C10,
		case when seqkey.plate_well = '3D07' then genotype end as 3D07
	from genotype join seqkey
		on genotype.dna_sample = seqkey.dna_id
		where call_id = 'GBS'
		and seqkey.sample like 'St37A'
		order by rs
);
drop view if exists tmp_gbsgenotype_pivot;
create view tmp_gbsgenotype_pivot as(
	select 
		rs, chrom,
		group_concat(3C10) as 3C10,
		group_concat(3D07) as 3D07
	from tmp_gbsgenotype
	group by chrom, rs
);
select 'rs', 'chrom', '3C10', '3D07'
UNION ALL
select * from tmp_gbsgenotype_pivot
into outfile '~/workdir/st37A.txt'
;
drop view if exists tmp_gbsgenotype;
create view tmp_gbsgenotype as(
	select rs, chrom, seqkey.sample, seqkey.plate_well, genotype,
		case when seqkey.plate_well = '3D10' then genotype end as 3D10,
		case when seqkey.plate_well = '3D08' then genotype end as 3D08
	from genotype join seqkey
		on genotype.dna_sample = seqkey.dna_id
		where call_id = 'GBS'
		and seqkey.sample like 'St38A'
		order by rs
);
drop view if exists tmp_gbsgenotype_pivot;
create view tmp_gbsgenotype_pivot as(
	select 
		rs, chrom,
		group_concat(3D10) as 3D10,
		group_concat(3D08) as 3D08
	from tmp_gbsgenotype
	group by chrom, rs
);
select 'rs', 'chrom', '3D10', '3D08'
UNION ALL
select * from tmp_gbsgenotype_pivot
into outfile '~/workdir/st38A.txt'
;
