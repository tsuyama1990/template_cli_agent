import argparse
import asyncio
from pathlib import Path
from ac_cdd_core.config import settings
from ac_cdd_core.services.git_ops import GitManager
from ac_cdd_core.services.aider_client import AiderClient
from ac_cdd_core.services.jules_client import JulesClient
from ac_cdd_core.sandbox import SandboxRunner

async def main():
    parser = argparse.ArgumentParser(description="Evaluate Code Review Capability")
    parser.add_argument("--pr", required=True, help="Pull Request URL")
    parser.add_argument("--session", required=True, help="Jules Session ID")
    args = parser.parse_args()
    
    git = GitManager()
    aider = AiderClient()
    jules = JulesClient()
    runner = SandboxRunner()
    
    print(f"=== Review Capability Test ===")
    print(f"PR: {args.pr}")
    print(f"Session: {args.session}")
    
    original_branch = None
    
    try:
        # 1. Checkout PR
        print("\n--- Checkout PR ---")
        try:
            original_branch = await git.get_current_branch()
        except:
            original_branch = "main"
            
        await git.checkout_pr(args.pr)
        current_branch = await git.get_current_branch()
        print(f"Checked out: {current_branch}")
        
        # 2. Identify Changed Files (Smart Audit)
        print("\n--- Identify Changes ---")
        # Assuming PR targets 'main'
        changed_files = await git.get_changed_files("main")
        print(f"Changed Files ({len(changed_files)}):")
        for f in changed_files:
            print(f" - {f}")
            
        # Context Files
        files_to_audit = []
        
        # Filter for code
        for f in changed_files:
            if f.endswith(".py") or f in ["ac_cdd_config.py", "pyproject.toml"]:
                if Path(f).exists():
                    files_to_audit.append(f)
                    
        # Add Architecture Doc (Context)
        arch_doc = Path(settings.paths.documents_dir) / "SYSTEM_ARCHITECTURE.md"
        if arch_doc.exists():
            print(f" - {arch_doc} (Context)")
            files_to_audit.append(str(arch_doc))
            
        if not files_to_audit:
            print("No relevant files to audit. Exiting.")
            return

        # 3. Conduct Review (Sanbox)
        print("\n--- Running Auditor (Aider) in Sandbox ---")
        
        # Load Auditor Instruction
        instr_path = Path(settings.paths.templates) / "AUDITOR_INSTRUCTION.md"
        if instr_path.exists():
            instruction = instr_path.read_text()
        else:
            instruction = "Conduct a strict code review."
            
        # Append context about Single Review
        instruction += "\n\n(Single Review Mode: Provide comprehensive feedback.)"
        
        output = await aider.run_audit(
            files=files_to_audit,
            instruction=instruction,
            runner=runner
        )
        
        print("\n--- Audit Result ---")
        # Print first 500 chars
        print(output[:500] + "..." if len(output) > 500 else output)
        
        # 4. Send API Feedback
        print(f"\n--- Sending Feedback to Session {args.session} ---")
        session_path = args.session
        if not session_path.startswith("sessions/"):
            session_path = f"sessions/{session_path}"
            
        # Use internal method for test
        await jules._send_message(session_path, output)
        print("Feedback sent successfully.")

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Cleanup
        print("\n--- Cleanup ---")
        try:
            if original_branch and original_branch != await git.get_current_branch():
                print(f"Restoring branch {original_branch}...")
                await git.runner.run_command(["git", "checkout", original_branch], check=False)
        except Exception as e:
            print(f"Cleanup failed: {e}")
            
        await runner.close()

if __name__ == "__main__":
    asyncio.run(main())
