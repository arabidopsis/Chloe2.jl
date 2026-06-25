
const rrn2product = Dict("rrn4.5"=>"4.5S rRNA","rrn5"=>"5S rRNA","rrn16"=>"16S rRNA","rrn23"=>"23S rRNA")

const rrn_model_lengths = Dict("rrn16_1" => 1491, "rrn23_1" => 2810, "rrn4.5_1" => 103, "rrn5_1" => 121)

function rrnsearch(targetfile::String; sensitivity = false)
    hmmpath = joinpath(chloe2models, "rrns", "all_rrns.hmm")
    results = PipeBuffer()
    nhmmer = which("nhmmer")
    cmd = sensitivity ? `$nhmmer --max -o /dev/null --tblout /dev/stdout --noali $hmmpath $targetfile` : `$nhmmer -o /dev/null --tblout /dev/stdout --noali $hmmpath $targetfile`
    run(cmd, devnull, results, stderr)
    return results
end

function parse_tbl(results::IOBuffer, glength::Integer)
    matches = FeatureMatch[]
    for line in readlines(results)
        startswith(line, "#") && continue
        bits = split(line, " ", keepempty=false)
        rrnstrand = bits[12][1]
        rrnstart = parse(Int, bits[7])
        rrnstart > glength && continue # don't retain matches starting in extension
        rrnstop = parse(Int, bits[8])
        model_from = parse(Int, bits[5])
        model_to = parse(Int, bits[6])
        fmstart = rrnstrand == '+' ? rrnstart : reverse_complement(rrnstart, glength)
        fmstart -= model_from - 1
        fm_stop = rrnstrand == '+' ? rrnstop : reverse_complement(rrnstop, glength)
        fmlength = fm_stop - fmstart + 1
        query = bits[3]
        fmlength += rrn_model_lengths[query] - model_to
        push!(matches, FeatureMatch(query, [query], rrnstrand, "rRNA", model_from, model_to,
            fmstart, fmlength, parse(Float64, bits[13])))
    end
    close(results)
    rationalise_matches!(matches, glength)
end

