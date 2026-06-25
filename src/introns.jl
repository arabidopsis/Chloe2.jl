const intron_model_lengths = Dict("atpF_2a" => 50,"atpF_2b" => 94,"clpP1_2a" => 103,"clpP1_2b" => 95,"clpP1_4a" => 170,"clpP1_4b" => 156,
"ndhA_2a" => 141,"ndhA_2b" => 105,"ndhB_2a" => 123,"ndhB_2b" => 158,"pafI_2a" => 184,"pafI_2b" => 95,
"pafI_4a" => 158,"pafI_4b" => 74,"petB_2a" => 184,"petB_2b" => 95,"petD_2a" => 210,"petD_2b" => 91,
"rpl16_2a" => 180,"rpl16_2b" => 94,"rpl2_2a" => 213,"rpl2_2b" => 123,"rpoC1_2a" => 152,"rpoC1_2b" => 93,
"rps12_2a" => 174,"rps12_2b" => 207,"rps12_3a" => 128,"rps12_3b" => 95,"rps12_5a" => 158,"rps12_5b" => 130,
"rps16_2a" => 162,"rps16_2b" => 74,"trnA-UGC_2a" => 137,"trnA-UGC_2b" => 138,"trnG-UCC_2a" => 156,
"trnG-UCC_2b" => 87,"trnI-GAU_2a" => 141,"trnI-GAU_2b" => 138,"trnK-UUU_2a" => 136,"trnK-UUU_2b" => 100,
"trnL-UAA_2a" => 107,"trnL-UAA_2b" => 77,"trnV-UAC_2a" => 219,"trnV-UAC_2b" => 172)

function cmsearch_intron(target::String, intron::String; sensitivity = false)  
    cmsearch = which("cmsearch")
    cmpath = joinpath(chloe2models, "introns", intron * "a.cm")
    resultsa = PipeBuffer()
    cmd = sensitivity ? `$cmsearch --max -o /dev/null --tblout /dev/stdout --toponly --notrunc --noali -Z 0.3 $cmpath $target` : `$cmsearch -o /dev/null --tblout /dev/stdout --toponly --hmmonly --notrunc --noali -Z 0.3 $cmpath $target`
    run(cmd, devnull, resultsa, stderr)

    cmpath = joinpath(chloe2models, "introns", intron * "b.cm")
    resultsb = PipeBuffer()
    cmd = sensitivity ? `$cmsearch --max -o /dev/null --tblout /dev/stdout --toponly --notrunc --noali -Z 0.3 $cmpath $target` : `$cmsearch -o /dev/null --tblout /dev/stdout --toponly --notrunc --noali -Z 0.3 $cmpath $target`
    run(cmd, devnull, resultsb, stderr)
    return resultsa, resultsb
end

function intronsearch(id::AbstractString, genome::CircularSequence, part, gene_model, tempfile::TempFile; sensitivity = false)
    #println(part)
    #println(gene_model)
    try
        leeway = 20 # arbitrary grey zone to account for slop in exon placement
        #defaults
        intron_strand = '+'
        donor_site = 1
        acceptor_site = length(genome)
        donor_idx = findlast(x -> only(partorder(x)) == part.order-1, gene_model)
        if ~isnothing(donor_idx)
            donor_match = gene_model[donor_idx]
            if donor_match.type ≠ "intron"    #to avoid searching for other half of trans-spliced intron
                donor_site = donor_match.target_from + donor_match.target_length - leeway
                intron_strand = donor_match.strand
            end
        end
        #println("$intron_strand\t$(string(donor_site))")
        acceptor_idx = findfirst(x -> only(partorder(x)) == part.order+1, gene_model)
        #println(acceptor_idx)
        if ~isnothing(acceptor_idx)
            acceptor_match = gene_model[acceptor_idx]
            if acceptor_match.type ≠ "intron"    #to avoid searching for other half of trans-spliced intron
                acceptor_site = acceptor_match.target_from + leeway
                intron_strand = acceptor_match.strand
            end
        end
        #println("$intron_strand\t$(string(acceptor_site))")
        if donor_site == 1 && acceptor_site < length(genome)
            donor_site = acceptor_site - 5000
        end
        if acceptor_site == length(genome) && donor_site > 1
            acceptor_site = donor_site + 5000
        end

        #don't search if search space is entire genome
        if donor_site == 1 && acceptor_site == length(genome); return missing; end

        search_seq = intron_strand == '+' ? genome[donor_site:acceptor_site] : reverse_complement(genome)[donor_site:acceptor_site]
        intron_name = "$(part.gene)_$(string(part.order))"
        fname = tempfilename(tempfile, "$id.$intron_name.$intron_strand.$(string(donor_site))-$(string(acceptor_site)).fa")
        #println("searching for $intron_name in $intron_strand.$(string(donor_site))-$(string(acceptor_site))")
        open(FASTA.Writer, fname) do writer
            write(writer, FASTA.Record("$id.$intron_name.$intron_strand.$(string(donor_site))-$(string(acceptor_site))", search_seq))
        end
        return cmsearch_intron(fname, "$intron_name"; sensitivity = sensitivity)
    catch
        println("failed intron search: $id $(part.gene)")
        return missing
    end
end

function parse_intron_tbl(results::IOBuffer, glength::Int)
    intron_matches = FeatureMatch[]
    hits = filter!(x -> ~startswith(x, "#"), readlines(results))
    isempty(hits) && return missing
    for hit in hits
        bits = split(hit, " ", keepempty=false)
        query = bits[3]
        m = match(r"(^[^.]+)(?:\.[0-9]+)?\.([^.]+)\.([+|-])\.(-?[0-9]+)-[0-9]+", bits[1]) # (?:\.[0-9]+)? is optional match to .1 version number on accession
        if isnothing(m)
            println("intronseqname: ", bits[1])
            println(hits)
        end
        id = m.captures[1]
        strand = first(m.captures[3])
        seqstart = parse(Int, m.captures[4])
        evalue = parse(Float64, bits[16])
        model_from = parse(Int, bits[6])
        model_to = parse(Int, bits[7])
        seq_from = parse(Int, bits[8])
        seq_to = parse(Int, bits[9])
        if query ∈ ["clpP1_4b", "rps12_5b"] # account for inset models
            model_to -= 20
            seq_to -= 20
        end
        target_from = seqstart + seq_from - model_from
        target_to = seqstart + seq_to - 1 + intron_model_lengths[query] - model_to
        target_length = target_to - target_from + 1
        # note that target coordinates are in stranded nucleotide coordinates
        push!(intron_matches, FeatureMatch(id, [query], strand, "intron", model_from, model_to, target_from, target_length, evalue))
    end
    close(results)
    @debug intronmatch
    rationalise_matches!(intron_matches, glength)
end

function parse_intron_tbls(results::Union{Missing, Tuple{IOBuffer, IOBuffer}}, glength)
    ismissing(results) && return missing
    matcha = parse_intron_tbl(results[1], glength)
    matchb = parse_intron_tbl(results[2], glength)
    if ismissing(matcha) || ismissing(matchb) || isempty(matcha) || isempty(matchb)
        return missing
    end
    rationalise_matches!(matcha, glength)
    rationalise_matches!(matchb, glength)
    merge_matches(first(matcha), first(matchb), glength)
end