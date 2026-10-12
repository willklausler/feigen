# Contributing

Bug reports, questions and pull requests are welcome.

- **Bugs and questions:** open an [issue](https://github.com/willklausler/feigen/issues) with the compiler, its version, the operating system, and the command line or a minimal program.
- **Pull requests:** fork, branch, and open a PR against `main`. CI must pass.

Before opening a PR:

1. Add tests for every change; run `fpm test --profile debug` and `fpm test --profile release`.
2. Run the examples and the app: `fpm run --example --all`, `fpm run`.
3. Document new or changed public procedures in `llms.txt`, `man/feigen.3` and `CHANGELOG.md`; `python3 tools/check_docs.py` must pass.

[AGENTS.md](AGENTS.md) describes the layout, conventions and compiler quirks in detail; it applies to human contributors too.
