
struct TempFile
    directory::String
    uuid::UUID
    ext::Vector{String}
    function TempFile(directory::String = ".")
        new(directory, uuid1(), [])
    end
end
function tempfilename(tempfile::TempFile, ext::String)
    push!(tempfile.ext, ext)
    joinpath(tempfile.directory, "$(tempfile.uuid).$(ext)")
end

function cleanfiles(tempfile::TempFile)
    while length(tempfile.ext) > 0
        ext = pop!(tempfile.ext)
        # don't add to tempfile.ext with tempfilename!!
        path = joinpath(tempfile.directory, "$(tempfile.uuid).$(ext)")
        rm(path; force = true)
    end
end

const EXES = Dict{String,String}()

function which(cmd::String)::String
    if cmd in keys(EXES)
        return EXES[cmd]
    end
    path = Sys.which(cmd)
    path === nothing && error("Required executable $cmd not found in PATH")
    EXES[cmd] = path
    return path
end

function missing_executables()::Vector{String}
    notfound = String[]
    for cmd in ["hmmsearch", "cmscan", "nhmmer", "cmsearch"]
        path = Sys.which(cmd)
        if path === nothing
            push!(notfound, cmd)
        else
            EXES[cmd] = path
        end
    end
    return notfound
end
function ensure_executables()
    notfound = missing_executables()
    if length(notfound) > 0
        s = length(notfound) == 1 ? " is" : "s are"
        p = length(notfound) == 1 ? "it" : "them"
        println(
            stderr,
            "The following required executable$(s) not in your PATH: \"$(join(notfound, ", "))\". " *
            "We can't continue without $(p). Please install $(p) and try again."
        )
        exit(0)
    end
end
