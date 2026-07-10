module Errors

export NotImplementedError
    
struct NotImplementedError <: Exception
    type::String
    NotImplementedError(type::String) = new(type)
    NotImplementedError(type) = new(string(typeof(type)))
end

function Base.showerror(io::IO, err::NotImplementedError)
    # print(io, "NotImplementedError: The method '$(err.func)' has not been implemented.")
    print(io, "NotImplementedError: No implementation for $(err.type).")
end

end # Errors