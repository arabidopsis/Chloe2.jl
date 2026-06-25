# unified struct to cover all matches of HMMs or CMs; target coordinates are from the 5' end of the strand.
mutable struct FeatureMatch
    target_id::String
    queryparts::Vector{AbstractString}
    strand::Char
    type::String
    model_from::Int
    model_to::Int
    target_from::Int
    target_length::Int
    evalue::Float64
end

function gene(fm::FeatureMatch)
    first(split(first(fm.queryparts), "_"))
end

function gene(gm::Vector{FeatureMatch})
    first(split(first(first(gm).queryparts), "_"))
end

function partorder(fm::FeatureMatch)
    sort!(fm.queryparts)
    parse.(Int, last.(split.(fm.queryparts, "_")))
end

Base.length(fm::FeatureMatch) = fm.target_length

function circularin(x::Integer, f::FeatureMatch, c::Integer)
    circularin(x, f.target_from, f.target_length, c)
end

function circularoverlap(f1::FeatureMatch, f2::FeatureMatch, c::Integer)
    if circularin(f1.target_from, f2, c)
        return circulardistance(f1.target_from, min(f1.target_from + f1.target_length -1, f2.target_from + f2.target_length - 1), c)
    elseif circularin(f2.target_from, f1, c)
        return circulardistance(f2.target_from, min(f1.target_from + f1.target_length -1, f2.target_from + f2.target_length - 1), c)
    end
    return 0
end

function circulardistance(m1::FeatureMatch, m2::FeatureMatch, c::Integer)
    return circulardistance(m1.target_from, m2.target_from, c)
end

function merge_matches(m1::FeatureMatch, m2::FeatureMatch, glength::Integer)
    sorted_features = sort([m1,m2]; by = x->(sort(x.queryparts)[1], x.model_from))
    m1 = sorted_features[1]
    m2 = sorted_features[2]
    # removes a and b suffixes from split intron models
    for i in eachindex(m1.queryparts)
        if ~isdigit(last(m1.queryparts[i]))
            m1.queryparts[i] = m1.queryparts[i][1:end-1]
        end
    end
    for i in eachindex(m2.queryparts)
        if ~isdigit(last(m2.queryparts[i]))
            m2.queryparts[i] = m2.queryparts[i][1:end-1]
        end
    end
    return FeatureMatch(m1.target_id, union(m1.queryparts, m2.queryparts), m1.strand, m1.type, m1.model_from, m2.model_to, m1.target_from,
        circulardistance(m1.target_from, m2.target_from + m2.target_length, glength), min(m1.evalue, m2.evalue))
end

function five2three(x, y, glength)
    circulardistance(y, x, glength) > circulardistance(x, y, glength)
end

# merge adjacent matches to same model
function rationalise_matches!(matches::Vector{FeatureMatch}, glength::Integer)::Vector{FeatureMatch}
    sort!(matches; by = x -> (x.strand, only(x.queryparts), x.evalue))
    todelete = Int[]
    bestmatchidx = 1
    for (i, m) in enumerate(matches)
        i == 1 && continue
        best = matches[bestmatchidx]
        if m.strand ≠ best.strand || only(m.queryparts) ≠ only(best.queryparts)
            bestmatchidx = i
            continue
        end
        #merge matches if they are hits to different parts of the same model
        modeloverlap = length(intersect(best.model_from:best.model_to, m.model_from:m.model_to)) / min(length(best.model_from:best.model_to), length(m.model_from:m.model_to))
        if genome_adjacent(best, m, glength) && (best.type ≠ "CDS" || in_frame(best, m)) && modeloverlap < 0.9
            matches[bestmatchidx] = merge_matches(best, m, glength)
            push!(todelete, i)
        elseif m.evalue > 1.05 * best.evalue #arbitrary 1.05 x threshold to cover variation in evalue even with identical targets
            push!(todelete, i)
        end
    end
    deleteat!(matches, todelete)
end

function group_duplicates(matches::Vector{FeatureMatch})
    query_dict = Dict{String, Vector{FeatureMatch}}()
    for item in matches
        query_value = item.query
        if haskey(query_dict, query_value)
            push!(query_dict[query_value], item)
        else
            query_dict[query_value] = [item]
        end
    end
    return query_dict
end

#Check if two CDS Features are in-frame
function in_frame(f1::FeatureMatch, f2::FeatureMatch)
    @assert f1.type == "CDS" && f2.type == "CDS"
    last(f1.target_id) == last(f2.target_id)
end

#should really check that frame is maintained when fixing junctions for CDS features
function fix_splice_junctions!(gene_model, glength)
    for (i,feature) in enumerate(gene_model)
        ismissing(feature) && continue
        if feature.type == "intron"
            if i-1 > 0
                previous_feature = gene_model[i - 1]
                if previous_feature.target_from + previous_feature.target_length ≠ feature.target_from # disagreement on boundary
                    if previous_feature.type == "tRNA"
                        feature.target_from = mod1(previous_feature.target_from + previous_feature.target_length, glength) # trust tRNA match over intron match
                    elseif previous_feature.type == "CDS"
                        previous_feature.target_length = circulardistance(previous_feature.target_from, feature.target_from, glength) # trust intron match over CDS match
                    end
                end
            end
            if i + 1 <= length(gene_model)
                next_feature = gene_model[i + 1]
                if feature.target_from + feature.target_length ≠ next_feature.target_from # disagreement on boundary
                    if next_feature.type == "tRNA"
                        feature.target_length = circulardistance(feature.target_from, next_feature.target_from, glength) # trust tRNA match over intron match
                    elseif next_feature.type == "CDS"
                        next_feature.target_from = mod1(feature.target_from + feature.target_length, glength) # trust intron match over CDS match
                    end
                end
            end
        end
    end
end

function genome_adjacent(m1::FeatureMatch, m2::FeatureMatch, glength)
    m1.strand ≠ m2.strand && return false
    #components within same transcription unit can be separated by over 3 kb (e.g. trnK exons), so use 4 kb as threshold for splitting transcription units
    min(closestdistance(m1.target_from + m1.target_length, m2.target_from, glength), closestdistance(m2.target_from + m2.target_length, m1.target_from, glength)) > 4000 && return false
    true
end

function addmatch2genemodel!(gene_model, feature_match, glength)
    for (i,exon) in enumerate(gene_model)
        if length(intersect(exon.queryparts, feature_match.queryparts)) > 0 # gene model already has a similar match
            if intersect(exon.model_from:exon.model_to, feature_match.model_from:feature_match.model_to) == 0
                gene_model[i] = merge_matches(exon, feature_match, glength)
                return true
            else
                return false
            end
        end
    end
    if isempty(gene_model) || genome_adjacent(last(gene_model), feature_match, glength)
        push!(gene_model, feature_match)
        return true
    end
    return false
end

function addmatches2genemodels!(part, gene_models, feature_matches, glength)
    matches = filter(x -> "$(part.gene)_$(part.order)" ∈ x.queryparts, feature_matches)
    sort!(matches; by = x -> (x.evalue, x.target_from))
    for m in matches
        placed = false
        for gm in gene_models
            placed = addmatch2genemodel!(gm, m, glength)
            placed && break
        end
        if ~placed
            push!(gene_models, [m]) # couldn't add to existing model, so start a new one
        end
    end
end

