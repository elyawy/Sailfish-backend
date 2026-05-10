# Sailfish — Emscripten build
# Prerequisites: source emsdk_env.sh before running make
#
# Output: sailfish.js + sailfish.wasm  (ES6 module)
# Usage:  make          — release build
#         make debug    — unoptimised + source maps
#         make clean

EMCC      := em++
OUTDIR    := dist
OUT       := $(OUTDIR)/sailfish.js

SRC       := wasm_bind.cpp

# Adjust paths if your layout differs
PHYLOLIB  := ../libs/emPhylolib/libPhylolib.a
INCLUDES  := -I../libs/pcg \
             -I../libs/Phylolib/includes

CXXFLAGS  := -std=c++17 \
             $(INCLUDES)

# -sMODULARIZE + EXPORT_ES6 → ES6 module  (import createSailfish from './sailfish.js')
# -sALLOW_MEMORY_GROWTH     → heap can grow; simulation sizes vary
# -lembind                  → Embind runtime
# -sEXPORT_NAME             → module factory function name in JS
EMFLAGS   := -lembind \
             -sMODULARIZE=1 \
             -sEXPORT_ES6=1 \
             -sEXPORT_NAME=createSailfish \
             -sALLOW_MEMORY_GROWTH=1 \
             -sENVIRONMENT=web

# ── Release ──────────────────────────────────────────────────────────────────
.PHONY: all
all: $(OUTDIR) $(OUT)

$(OUT): $(SRC) $(PHYLOLIB)
	$(EMCC) $(CXXFLAGS) -O2 $(SRC) $(PHYLOLIB) $(EMFLAGS) -o $(OUT)

# ── Debug (no optimisation, DWARF source maps) ────────────────────────────────
.PHONY: debug
debug: $(OUTDIR)
	$(EMCC) $(CXXFLAGS) -O0 -g $(SRC) $(PHYLOLIB) $(EMFLAGS) \
	    -sASSERTIONS=2 \
	    -o $(OUTDIR)/sailfish.debug.js

# ── Helpers ───────────────────────────────────────────────────────────────────
$(OUTDIR):
	mkdir -p $(OUTDIR)

.PHONY: clean
clean:
	rm -rf $(OUTDIR)
