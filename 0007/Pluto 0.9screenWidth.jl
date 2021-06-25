### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 3ecc9360-d570-11eb-0a04-ad57130a3b16
using HypertextLiteral

# ╔═╡ a6033a20-bdb6-41bd-bb31-41e827f756e0
@bind screenWidth @htl("""
	<div>
	<script>
		var div = currentScript.parentElement
		div.value = screen.width
	</script>
	</div>
""")

# ╔═╡ c971f8aa-833f-4ec4-99a6-c3e10437c57a
begin
	cellWidth= min(1000, screenWidth*0.9)
	@htl("""
		<style>
			pluto-notebook {
				margin: auto;
				width: $(cellWidth)px;
			}
		</style>
	""")
end

# ╔═╡ Cell order:
# ╠═3ecc9360-d570-11eb-0a04-ad57130a3b16
# ╠═a6033a20-bdb6-41bd-bb31-41e827f756e0
# ╠═c971f8aa-833f-4ec4-99a6-c3e10437c57a
