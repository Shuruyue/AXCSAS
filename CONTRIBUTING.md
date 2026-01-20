# Contributing to AXCSAS

Thank you for your interest in contributing to AXCSAS! This document provides guidelines and instructions for contributors.

## Table of Contents

1. [Code of Conduct](#code-of-conduct)
2. [Getting Started](#getting-started)
3. [Development Setup](#development-setup)
4. [Making Changes](#making-changes)
5. [Testing](#testing)
6. [Code Style](#code-style)
7. [Submitting Changes](#submitting-changes)
8. [Adding New Features](#adding-new-features)

---

## Code of Conduct

Be respectful, constructive, and professional in all interactions.

---

## Getting Started

**Prerequisites**:
- Python ≥ 3.8
- Git
- Basic knowledge of XRD analysis

**Fork and Clone**:
```bash
# Fork the repository on GitHub, then:
git clone https://github.com/YOUR_USERNAME/AXCSAS.git
cd AXCSAS
```

---

## Development Setup

**1. Create virtual environment**:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

**2. Install in editable mode**:
```bash
pip install -e ".[dev]"
```

**3. Install development tools**:
```bash
pip install pytest ruff black mypy
```

**4. Verify installation**:
```bash
pytest tests/
```

---

## Making Changes

**1. Create a branch**:
```bash
git checkout -b feature/your-feature-name
# or
git checkout -b bugfix/issue-number
```

**2. Make your changes following our code standards**

**3. Write or update tests**

**4. Run the test suite**:
```bash
pytest tests/ -v
```

---

## Testing

**Run all tests**:
```bash
pytest tests/
```

**Run specific test**:
```bash
pytest tests/test_scherrer.py -v
```

**Check coverage**:
```bash
pytest --cov=axcsas --cov-report=html
```

**Writing Tests**:
- Place tests in `tests/` directory
- Use descriptive test names: `test_scherrer_calculates_size_correctly()`
- Include edge cases and error conditions
- Target >90% coverage for new code

---

## Code Style

**Python Style**:
- Follow PEP 8
- Use type hints for function signatures
- Maximum line length: 100 characters
- Use double quotes for strings

**Format code automatically**:
```bash
black axcsas/
```

**Check linting**:
```bash
ruff check axcsas/
```

**Type checking**:
```bash
mypy axcsas/
```

**Documentation**:
- All public functions must have docstrings
- Use Google-style docstrings
- Include Args, Returns, and Examples sections
- Bilingual comments (English + Chinese) where helpful

**Example**:
```python
def calculate_size(
    two_theta: float,
    fwhm: float,
    wavelength: float = CU_KA1
) -> float:
    """
    Calculate crystallite size using Scherrer equation.
    
    Args:
        two_theta: Diffraction angle in degrees
        fwhm: Full width at half maximum in degrees
        wavelength: X-ray wavelength in Angstroms
        
    Returns:
        Crystallite size in nanometers
        
    Example:
        >>> size = calculate_size(43.3, 0.35)
        >>> print(f"{size:.1f} nm")
    """
    ...
```

---

## Submitting Changes

**1. Commit your changes**:
```bash
git add .
git commit -m "feat: Add support for aluminum analysis"
```

**Commit Message Format**:
```
<type>: <subject>

<body>
```

**Types**:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting)
- `refactor`: Code refactoring
- `test`: Adding/updating tests
- `chore`: Maintenance tasks

**2. Push to your fork**:
```bash
git push origin feature/your-feature-name
```

**3. Create Pull Request**:
- Go to the original repository on GitHub
- Click "New Pull Request"
- Select your branch
- Fill in the PR template
- Wait for review

**PR Checklist**:
- [ ] Tests pass locally
- [ ] Code follows style guidelines
- [ ] Documentation updated
- [ ] CHANGELOG.md updated
- [ ] No merge conflicts

---

## Adding New Features

**Adding New Materials**:

1. Add constants to `core/constants.py` or create new module like `core/aluminum_crystal.py`
2. Update analysis pipeline to support new material
3. Add tests in `tests/test_new_material.py`
4. Document in API reference

**Adding New Analysis Methods**:

1. Create new module in `axcsas/methods/`
2. Follow existing patterns (dataclasses for results, analyzer classes)
3. Add comprehensive docstrings
4. Include academic citations
5. Add visualization support
6. Write unit tests
7. Update CLI if needed

**Example Structure**:
```python
# axcsas/methods/new_method.py

from dataclasses import dataclass
from typing import Optional

@dataclass
class NewMethodResult:
    """Result of new analysis method."""
    value: float
    uncertainty: float
    is_valid: bool

class NewMethodAnalyzer:
    """
    Analyzer for new XRD analysis method.
    
    Reference:
        Author, A. (Year). Title. Journal, vol(issue), pages.
    """
    
    def analyze(self, data) -> NewMethodResult:
        """Run analysis."""
        ...
```

---

## Questions?

- Open an issue for bugs or feature requests
- Start a discussion for questions
- Email: [maintainer email]

---

**Thank you for contributing to AXCSAS!** 🎉
