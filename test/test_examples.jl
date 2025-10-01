module Examples
    
examples_main = []

TEST_EX = true
for (root, dirs, files) in walkdir(joinpath(dirname(@__FILE__), "../examples"))
    for file in files
        include(joinpath(root, file))
        push!(examples_main, (name=file, main=Example.main))
    end
end
TEST_EX = false

function test()
    for ex_item in examples_main
        try
            ex_item.main(true)
        catch e
            println("Error in $(ex_item.name): ", e)
            println("Stacktrace:")
            showerror(stdout, e, catch_backtrace())
            return false
        end
    end

    return true
end

end # Examples