# Keratia Codespaces Cheatsheet

This repo is intended to run in a GitHub Codespace using the default Python environment, without creating a separate virtual environment.

## 1. Install or update Python requirements

From the repository root:

```bash
python -m pip install --user -r requirements.txt
```

If Jupyter does not immediately see the same Python environment, register a named kernel:

```bash
python -m ipykernel install --user --name keratia --display-name "Python (keratia)"
```

Then select `Python (keratia)` as the notebook kernel.

## 2. Check where you are

```bash
pwd
ls
```

You should see files such as:

```text
requirements.txt
cheatsheet.md
```

## 3. Check Git status

```bash
git status
```

This tells you whether your Codespace has local changes that are not committed.

## 4. Fetch updates from the remote

Terminology note: `origin` is the default name of the remote repository. `fetch` downloads new commits and branch information, but does not update your working files by itself.

```bash
git fetch origin
```

To see whether your local branch is behind the remote branch:

```bash
git status
```

## 5. Pull updates from the remote branch

If you want your local branch to actually update to match the remote branch, use `pull`.

For the main branch:

```bash
git pull origin main
```

In this workflow, if ChatGPT commits directly to the GitHub repo, you should usually run:

```bash
git pull origin main
```

inside your Codespace before testing the latest changes.

## 6. If you have local edits before pulling

Check them first:

```bash
git status
```

If you want to keep your edits but temporarily move them aside:

```bash
git stash push -m "temporary local edits"
git pull origin main
git stash pop
```

If you want to discard local changes completely, be careful:

```bash
git restore .
git pull origin main
```

## 7. Useful notebook command

Launch Jupyter from the repo root if needed:

```bash
jupyter notebook
```

In GitHub Codespaces, VS Code's built-in notebook interface is usually easier than manually launching a notebook server.
