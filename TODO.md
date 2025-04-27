# TODO

## 🚀 Upcoming Features
- [ ] Expose a high-level `Hamiltonian` builder for spin systems

## 📚 Documentation
- [ ] Expand docstrings with usage examples for each function
- [ ] Generate and publish Sphinx HTML docs (enable ReadTheDocs webhook)
- [ ] Add more example notebooks under `examples/`
- [ ] Create a usage diagram in the README

## ✅ Testing
- [ ] Write pytest tests for:
  - `tlm`, `unit_tlm`, `number_of_states`
  - `getIy`, `getIz`, `getIp`, `getIm`, `getTlm`, `getEf`, `getIxf`, `getIyf`, `getIzf`, `getIpf`, `getImf`
  - utility functions `mypow`, `fac`
- [ ] Add a CI job to install and test built wheels in a fresh virtualenv
- [ ] Integrate coverage reporting and add badge to the README

## 🔧 CI & Automation
- [ ] Finalize pre-commit hooks (`black`, `flake8`, `isort`, `mypy`)
- [ ] Configure GitHub Actions `test_wheels` job to verify wheel install
- [ ] Add cache for Cython build artifacts (if beneficial)

## 🧹 Housekeeping
- [ ] Update `CHANGELOG.md` for the next release
- [ ] Bump version in `pyproject.toml` or `setup.py` before tagging
- [ ] Review and consolidate any TODO comments in code
