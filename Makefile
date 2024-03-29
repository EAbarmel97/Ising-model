include .env

JULIA_DEPOT_PATH := $(shell pwd)/.julenv

define julia_env
	julia --project=$(JULIA_DEPOT_PATH) $(1)
endef

define icn_julia_env
	$(ICN_JULIA_BIN) --project=$(ICN_JULIA_DEPOT_PATH) $(1)
endef

createprojectdir:
	@mkdir -p $(JULIA_DEPOT_PATH)

juliaenv:
	@$(call julia_env,$(ARGS))
	
instantiatenv: createprojectdir
	@cp Project.toml $(JULIA_DEPOT_PATH)/Project.toml
	@$(call julia_env,-e 'using Pkg; Pkg.resolve(); Pkg.instantiate()')

icncreateprojectdir:
	@mkdir -p $(ICN_JULIA_DEPOT_PATH)	

icnjuliaenv:
	@$(call icn_julia_env,$(ARGS))

icninstantiatenv: icncreateprojectdir
	@cp Project.toml $(ICN_JULIA_DEPOT_PATH)/Project.toml
	@$(call icn_julia_env,-e 'using Pkg; Pkg.resolve(); Pkg.instantiate()')

.PHONY: createprojectdir juliaenv instantiatenv icncreateprojectdir icnjuliaenv icninstantiatenv